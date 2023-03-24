using Oceananigans
using Oceananigans.TimeSteppers: tick!
using GLMakie

arch = CPU()
oNx = oNy = 256
aNx = aNy = 256
oNz = 1
aNz = 1

ocean_Δt = 0.01
atmos_Δt = 0.05
coupling_Δt = 0.05

ρ_ocean = 10.0
ρ_atmos = 1.0
Cᴰ = 1.0

ocean_grid = RectilinearGrid(arch, size=(oNx, oNy, oNz), x=(0, 2π), y=(0, 2π), z=(0, 1))
atmos_grid = RectilinearGrid(arch, size=(aNx, aNy, aNz), x=(0, 2π), y=(0, 2π), z=(0, 1))

# Regrid to atmospheric grid for exchanges
surface_grid = RectilinearGrid(arch, size=(aNx, aNy, 1), x=(0, 2π), y=(0, 2π), z=(0, 1))

ocean_τx = XFaceField(ocean_grid)
atmos_τx = XFaceField(atmos_grid)
ocean_τy = YFaceField(ocean_grid)
atmos_τy = YFaceField(atmos_grid)

ocean_u_bcs = FieldBoundaryConditions(top    = FluxBoundaryCondition(interior(ocean_τx, :, :, 1)))
ocean_v_bcs = FieldBoundaryConditions(top    = FluxBoundaryCondition(interior(ocean_τy, :, :, 1)))
atmos_u_bcs = FieldBoundaryConditions(bottom = FluxBoundaryCondition(interior(atmos_τx, :, :, 1)))
atmos_v_bcs = FieldBoundaryConditions(bottom = FluxBoundaryCondition(interior(atmos_τy, :, :, 1)))

atmos_model = NonhydrostaticModel(; grid=atmos_grid, advection=WENO(),
                                  boundary_conditions = (u=atmos_u_bcs, v=atmos_v_bcs))  
                    
ocean_model = NonhydrostaticModel(; grid=ocean_grid, advection=WENO(),
                                  boundary_conditions = (u=ocean_u_bcs, v=ocean_v_bcs))  

ϵ(x, y, z) = 2rand() - 1
set!(ocean_model, u=ϵ, v=ϵ)

atmos_simulation = Simulation(atmos_model; Δt=atmos_Δt, stop_time=0)
ocean_simulation = Simulation(ocean_model; Δt=ocean_Δt, stop_time=0)

atmos_wizard = TimeStepWizard(cfl=0.2, max_change=1.1)
atmos_simulation.callbacks[:wizard] = Callback(atmos_wizard, IterationInterval(10))

ocean_wizard = TimeStepWizard(cfl=0.2, max_change=1.1)
ocean_simulation.callbacks[:wizard] = Callback(ocean_wizard, IterationInterval(10))

ocean_simulation.stop_time = 1
run!(ocean_simulation)
ocean_model.clock.time = 0.0
ocean_model.clock.iteration = 0

struct FluxExchanger{O, A, OF, AF}
    ocean_state :: O
    atmos_state :: A
    ocean_surface_fluxes:: OF
    atmos_surface_fluxes:: AF
    drag_coefficient :: Float64
    ocean_density :: Float64
    atmos_density :: Float64
end

struct CouplerState{U, V}
    u :: U
    v :: V
end

struct SurfaceFluxes{X, Y}
    τx :: X
    τy :: Y
end

ocean_surface_fluxes   = SurfaceFluxes(ocean_τx, ocean_τy)
atmos_surface_fluxes   = SurfaceFluxes(atmos_τx, atmos_τy)

coupler_ocean_surface_u = XFaceField(surface_grid)
coupler_ocean_surface_v = YFaceField(surface_grid)
coupler_atmos_surface_u = XFaceField(surface_grid)
coupler_atmos_surface_v = YFaceField(surface_grid)

ocean_state = CouplerState(coupler_ocean_surface_u, coupler_ocean_surface_v)
atmos_state = CouplerState(coupler_atmos_surface_u, coupler_atmos_surface_v)

flux_exchanger = FluxExchanger(ocean_state,
                               atmos_state,
                               ocean_surface_fluxes,
                               atmos_surface_fluxes,
                               Cᴰ, ρ_ocean, ρ_atmos)

struct CoupledSimulation{F, O, A}
    flux_exchanger :: F
    ocean_simulation :: O
    atmos_simulation :: A
    clock :: Clock{Float64}
end

clock = Clock(time=0.0)
cs = CoupledSimulation(flux_exchanger, ocean_simulation, atmos_simulation, clock)

function time_step!(cs::CoupledSimulation, Δt)

    # Update the coupler state
    oNx, oNy, oNz = size(cs.ocean_simulation.model.grid)
    ocean_surface_u = interior(cs.ocean_simulation.model.velocities.u, :, :, oNz)
    ocean_surface_v = interior(cs.ocean_simulation.model.velocities.v, :, :, oNz)

    aNx, aNy, aNz = size(cs.atmos_simulation.model.grid)
    atmos_surface_u = interior(cs.atmos_simulation.model.velocities.u, :, :, 1)
    atmos_surface_v = interior(cs.atmos_simulation.model.velocities.v, :, :, 1)

    fe = cs.flux_exchanger
    set!(fe.ocean_state.u, ocean_surface_u)
    set!(fe.ocean_state.v, ocean_surface_v)

    set!(fe.atmos_state.u, atmos_surface_u)
    set!(fe.atmos_state.v, atmos_surface_v)

    # Compute surface fluxes
    uo = interior(fe.ocean_state.u)
    vo = interior(fe.ocean_state.v)
    ua = interior(fe.atmos_state.u)
    va = interior(fe.atmos_state.v)
    τax = interior(fe.atmos_surface_fluxes.τx)
    τay = interior(fe.atmos_surface_fluxes.τy)
    τox = interior(fe.ocean_surface_fluxes.τx)
    τoy = interior(fe.ocean_surface_fluxes.τy)

    ρa = fe.atmos_density
    ρo = fe.ocean_density
    Cᴰ = fe.drag_coefficient

    @. τax = Cᴰ * (uo - ua)^2
    @. τay = Cᴰ * (vo - va)^2 

    # Update the surface fluxes in component models
    @. τox = ρa / ρo * τax
    @. τoy = ρa / ρo * τay

    # Run component models
    cs.ocean_simulation.stop_time += Δt
    run!(cs.ocean_simulation)

    cs.atmos_simulation.stop_time += Δt
    run!(cs.atmos_simulation)

    tick!(cs.clock, Δt)

    return nothing
end

ocean_simulation.verbose = false
atmos_simulation.verbose = false

uo, vo, wo = ocean_model.velocities
ua, va, wa = atmos_model.velocities
ζo = compute!(Field(∂x(vo) - ∂y(uo)))
ζa = compute!(Field(∂x(va) - ∂y(ua)))

ζot = []
ζat = []

push!(ζot, Array(interior(compute!(ζo), :, :, 1)))
push!(ζat, Array(interior(compute!(ζa), :, :, 1)))

t = [0.0]
tf = 5.0

while t[end] < tf
    time_step!(cs, atmos_Δt)
    push!(t, t[end] + atmos_Δt)
    if n % 10 == 0
        push!(ζot, Array(interior(compute!(ζo), :, :, 1)))
        push!(ζat, Array(interior(compute!(ζa), :, :, 1)))
        @info string("t: ", t[n],
                     ", atmos iteration: ", iteration(atmos_simulation),
                     ", atmos Δt : ", atmos_simulation.Δt,
                     ", ocean iteration: ", iteration(ocean_simulation),
                     ", ocean Δt: ", ocean_simulation.Δt)
    end
end

fig = Figure()

axo = Axis(fig[1, 1], title="Ocean")
axa = Axis(fig[1, 2], title="Atmosphere")
slider = Slider(fig[2, 1:2], range=1:Nt, startvalue=1)
n = slider.value

ζon = @lift ζot[$n]
ζan = @lift ζat[$n]

ζolim = 10
ζalim = 0.5

heatmap!(axo, ζon, colormap=:balance, colorrange=(-ζolim, ζolim))
heatmap!(axa, ζan, colormap=:balance, colorrange=(-ζalim, ζalim))

#display(fig)

