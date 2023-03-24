using Oceananigans
using Oceananigans.Units
using DoubleTrouble: CoupledAtmosphereOceanModel, relative_vertical_vorticity
using GLMakie

# Defint the grid
arch = CPU()
oNx = aNx = 360
oNy = aNy = 120
oNz = 1
aNz = 1
latitude = (-60, 60)
longitude = (0, 360)
halo = (5, 5, 5)
topology = (Periodic, Bounded, Bounded)
z = (0, Lz)
Lz = 1000

# Characteristic initial velocities
U_ocean = 1.0
U_atmos = 20.0

# Simulation parameters
ocean_Δt = 20minutes
atmos_Δt = 10minutes
coupling_Δt = 2 * ocean_Δt
stop_time = 1day
output_frequency = 2hours 

# Simple setup for atmosphere and ocean
model_setup(grid) = (tracers = nothing,
                     buoyancy = nothing,
                     coriolis = HydrostaticSphericalCoriolis(),
                     momentum_advection = VectorInvariant(vorticity_scheme   = WENO(grid),
                                                          divergence_scheme  = WENO(grid),
                                                          vertical_scheme    = WENO(grid)))

# Atmospheric model
atmos_grid = LatitudeLongitudeGrid(arch, size=(aNx, aNy, aNz); latitude, longitude, z, halo, topology)
atmos_τx = XFaceField(atmos_grid, indices=(:, :, 1))
atmos_τy = YFaceField(atmos_grid, indices=(:, :, 1))
atmos_u_bcs = FieldBoundaryConditions(bottom = FluxBoundaryCondition(interior(atmos_τx, :, :, 1)))
atmos_v_bcs = FieldBoundaryConditions(bottom = FluxBoundaryCondition(interior(atmos_τy, :, :, 1)))

atmos_model = HydrostaticFreeSurfaceModel(; grid = atmos_grid,
                                          model_setup(atmos_grid)...,
                                          boundary_conditions = (u=atmos_u_bcs, v=atmos_v_bcs))  

ϵ(x, y, z) = U_atmos * (2rand() - 1)
set!(atmos_model, u=ϵ, v=ϵ)

atmos_simulation = Simulation(atmos_model; Δt=atmos_Δt)
atmos_outputs = (; ζ = relative_vertical_vorticity(atmos_model))
atmos_simulation.output_writers[:jld2] = JLD2OutputWriter(atmos_model, atmos_outputs,
                                                          schedule = TimeInterval(output_frequency),
                                                          filename = "atmos_output.jld2",
                                                          overwrite_existing = true)

# Ocean model
ocean_grid = LatitudeLongitudeGrid(arch, size=(oNx, oNy, oNz); latitude, longitude, z, halo, topology)
ocean_τx = XFaceField(ocean_grid, indices=(:, :, oNz))
ocean_τy = YFaceField(ocean_grid, indices=(:, :, oNz))
ocean_u_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(interior(ocean_τx, :, :, 1)))
ocean_v_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(interior(ocean_τy, :, :, 1)))
ocean_model = HydrostaticFreeSurfaceModel(; grid = ocean_grid,
                                  model_setup(atmos_grid)...,
                                  boundary_conditions = (u=ocean_u_bcs, v=ocean_v_bcs))  


ϵ(x, y, z) = U_ocean * (2rand() - 1)
set!(ocean_model, u=ϵ, v=ϵ)
ocean_simulation = Simulation(ocean_model; Δt=ocean_Δt)
ocean_outputs = (; ζ = relative_vertical_vorticity(ocean_model))
ocean_simulation.output_writers[:jld2] = JLD2OutputWriter(ocean_model, ocean_outputs,
                                                          schedule = TimeInterval(output_frequency),
                                                          filename = "ocean_output.jld2",
                                                          overwrite_existing = true)

# Coupled model
surface_grid = LatitudeLongitudeGrid(arch, size=(aNx, aNy); latitude, longitude, topology=(Periodic, Bounded, Flat))

coupled_model = CoupledAtmosphereOceanModel(surface_grid,
                                            ocean = ocean_simulation,
                                            atmosphere = atmos_simulation)

coupled_simulation = Simulation(coupled_model; Δt=coupling_Δt, stop_time)

progress(sim) = @info string("t: ", prettytime(sim),
                             " coupled iter: ", iteration(sim),
                             ", atmos iter: ", iteration(sim.model.components.atmos),
                             ", atmos Δt : ", prettytime(sim.model.components.atmos.Δt),
                             ", ocean iter: ", iteration(sim.model.components.ocean),
                             ", ocean Δt : ", prettytime(sim.model.components.ocean.Δt))
 
coupled_simulation.callbacks[:progress] = Callback(progress, IterationInterval(10))

run!(coupled_simulation)

ζot = FieldTimeSeries("ocean_output.jld2", "ζ")
ζat = FieldTimeSeries("atmos_output.jld2", "ζ")
t = ζot.times
Nt = length(t)

fig = Figure(resolution=(1600, 800))

axo = Axis(fig[1, 1], xlabel="x (km)", ylabel="y (km)", title="Ocean vorticity")
axa = Axis(fig[1, 2], xlabel="x (km)", ylabel="y (km)", title="Atmosphere vorticity")
slider = Slider(fig[2, 1:2], range=1:Nt, startvalue=1)
n = slider.value

title = @lift string("Two-dimensional coupled atmosphere-ocean at t = ",
                     prettytime(t[$n]))
Label(fig[0, 1:2], title)

ζon = @lift interior(ζot[$n], :, :, 1)
ζan = @lift interior(ζat[$n], :, :, 1)

ζolim = maximum(abs, ζot[1]) / 10
ζalim = maximum(abs, ζat[1]) / 10

x, y, z = nodes(ζot)

x = x ./ 1e3
y = y ./ 1e3

heatmap!(axo, x, y, ζon, colormap=:balance, colorrange=(-ζolim, ζolim))
heatmap!(axa, x, y, ζan, colormap=:balance, colorrange=(-ζalim, ζalim))

display(fig)

record(fig, "two_dimensional_atmosphere_ocean.mp4", 1:Nt) do nn
    n[] = nn
end

