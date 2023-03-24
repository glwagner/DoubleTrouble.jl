using Oceananigans
using Oceananigans.Units
using DoubleTrouble: CoupledAtmosphereOceanModel, relative_vertical_vorticity
using GLMakie

arch = CPU()

# Resolution
oNx = oNy = 128
aNx = aNy = 128
oNz = 1
aNz = 1

# Domain
Lx = Ly = 1000kilometers
Lz = 1kilometer

# Characteristic initial velocities
U_ocean = 2.0
U_atmos = 10.0

ocean_Δt = 0.1 * (Lx / oNx) / U_ocean
atmos_Δt = 0.1 * (Lx / aNx) / U_atmos
coupling_Δt = 2 * ocean_Δt

# Atmospheric model
atmos_grid = RectilinearGrid(arch, size=(aNx, aNy, aNz), x=(0, Lx), y=(0, Ly), z=(0, Lz))
atmos_τx = XFaceField(atmos_grid, indices=(:, :, 1))
atmos_τy = YFaceField(atmos_grid, indices=(:, :, 1))
atmos_u_bcs = FieldBoundaryConditions(bottom = FluxBoundaryCondition(interior(atmos_τx, :, :, 1)))
atmos_v_bcs = FieldBoundaryConditions(bottom = FluxBoundaryCondition(interior(atmos_τy, :, :, 1)))

atmos_model = NonhydrostaticModel(; grid = atmos_grid,
                                  advection = WENO(),
                                  boundary_conditions = (u=atmos_u_bcs, v=atmos_v_bcs))  

ϵ(x, y, z) = U_atmos * (2rand() - 1)
set!(atmos_model, u=ϵ, v=ϵ)

atmos_simulation = Simulation(atmos_model; Δt=atmos_Δt)
atmos_wizard = TimeStepWizard(cfl=0.2, max_change=1.1)
atmos_simulation.callbacks[:wizard] = Callback(atmos_wizard, IterationInterval(10))
atmos_outputs = (; ζ = relative_vertical_vorticity(atmos_model))
atmos_simulation.output_writers[:jld2] = JLD2OutputWriter(atmos_model, atmos_outputs,
                                                          schedule = IterationInterval(10),
                                                          filename = "atmos_output.jld2",
                                                          overwrite_existing = true)

# Ocean model
ocean_grid = RectilinearGrid(arch, size=(oNx, oNy, oNz), x=(0, Lx), y=(0, Ly), z=(0, Lz))
ocean_τx = XFaceField(ocean_grid, indices=(:, :, oNz))
ocean_τy = YFaceField(ocean_grid, indices=(:, :, oNz))
ocean_u_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(interior(ocean_τx, :, :, 1)))
ocean_v_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(interior(ocean_τy, :, :, 1)))
ocean_model = NonhydrostaticModel(; grid=ocean_grid, advection=WENO(),
                                  boundary_conditions = (u=ocean_u_bcs, v=ocean_v_bcs))  


ϵ(x, y, z) = U_ocean * (2rand() - 1)
set!(ocean_model, u=ϵ, v=ϵ)
ocean_simulation = Simulation(ocean_model; Δt=ocean_Δt)
ocean_wizard = TimeStepWizard(cfl=0.2, max_change=1.1)
ocean_simulation.callbacks[:wizard] = Callback(ocean_wizard, IterationInterval(10))
ocean_outputs = (; ζ = relative_vertical_vorticity(ocean_model))
ocean_simulation.output_writers[:jld2] = JLD2OutputWriter(ocean_model, ocean_outputs,
                                                          schedule = IterationInterval(10),
                                                          filename = "ocean_output.jld2",
                                                          overwrite_existing = true)

# Coupled model
surface_grid = RectilinearGrid(arch, size=(aNx, aNy), x=(0, Lx), y=(0, Ly),
                               topology=(Periodic, Periodic, Flat))

coupled_model = CoupledAtmosphereOceanModel(surface_grid,
                                            ocean = ocean_simulation,
                                            atmosphere = atmos_simulation)

coupled_simulation = Simulation(coupled_model, Δt=coupling_Δt, stop_iteration=500)

progress(sim) = @info string("t: ", prettytime(sim),
                             ", atmos iteration: ", iteration(sim.model.components.atmos),
                             ", atmos Δt : ",       prettytime(sim.model.components.atmos.Δt),
                             ", ocean iteration: ", iteration(sim.model.components.ocean),
                             ", ocean Δt : ",       prettytime(sim.model.components.ocean.Δt))
 
coupled_simulation.callbacks[:progress] = Callback(progress, IterationInterval(10))

run!(coupled_simulation)

ζot = FieldTimeSeries("ocean_output.jld2", "ζ")
ζat = FieldTimeSeries("atmos_output.jld2", "ζ")
t = ζot.times
Nt = length(t)

fig = Figure()

axo = Axis(fig[1, 1], title="Ocean")
axa = Axis(fig[1, 2], title="Atmosphere")
slider = Slider(fig[2, 1:2], range=1:Nt, startvalue=1)
n = slider.value

ζon = @lift interior(ζot[$n], :, :, 1)
ζan = @lift interior(ζat[$n], :, :, 1)

ζolim = maximum(abs, ζot[1]) / 10
ζalim = maximum(abs, ζat[1]) / 10

heatmap!(axo, ζon, colormap=:balance, colorrange=(-ζolim, ζolim))
heatmap!(axa, ζan, colormap=:balance, colorrange=(-ζalim, ζalim))

display(fig)

