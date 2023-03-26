using DoubleTrouble: CoupledAtmosphereOceanModel
using DoubleTrouble.Utils: plane_simulation
using DoubleTrouble.Fluxes: Flux, AirSeaQuadraticDrag, AirSeaTracerTransfer

using Oceananigans
using Oceananigans.Units
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
U_ocean = 1.0
U_atmos = 20.0

ocean_Δt = 0.1 * (Lx / oNx) / U_ocean
atmos_Δt = 0.1 * (Lx / aNx) / U_atmos
coupling_Δt = 2 * ocean_Δt

x = (0, Lx)
y = (0, Ly)
tracers = (:T, :CO₂)
atmos = plane_simulation(size=(aNx, aNy); x, y, tracers, filename="atmos.jld2", interface=:bottom)
ocean = plane_simulation(size=(oNx, oNy); x, y, tracers, filename="ocean.jld2", interface=:top)

ϵ(x, y, z) = U_atmos * (2rand() - 1)
set!(atmos.model, u=ϵ, v=ϵ)

ϵ(x, y, z) = U_ocean * (2rand() - 1)
set!(ocean.model, u=ϵ, v=ϵ)

@show quadratic_drag = AirSeaQuadraticDrag()
@show heat_transfer = AirSeaTracerTransfer(:T, 1e-2)
@show co2_transfer = AirSeaTracerTransfer(:CO₂, 1e-3)

grid = atmos.model.grid
momentum_flux = Flux(quadratic_drag, grid)
heat_flux = Flux(heat_transfer, grid)
carbon_flux = Flux(co2_transfer, grid)

model = CoupledAtmosphereOceanModel(ocean = ocean,
                                    atmosphere = atmos,
                                    fluxes = (momentum_flux, heat_flux, carbon_flux))

#=
coupled_simulation = Simulation(coupled_model, Δt=coupling_Δt, stop_iteration=10)

progress(sim) = @info string("t: ", prettytime(sim),
                             " coupled iter: ", iteration(sim),
                             ", atmos iter: ", iteration(sim.model.components.atmos),
                             ", atmos Δt : ",       prettytime(sim.model.components.atmos.Δt),
                             ", ocean iter: ", iteration(sim.model.components.ocean),
                             ", ocean Δt : ",       prettytime(sim.model.components.ocean.Δt))
 
coupled_simulation.callbacks[:progress] = Callback(progress, IterationInterval(10))

run!(coupled_simulation)

ζot = FieldTimeSeries("ocean.jld2", "ζ")
ζat = FieldTimeSeries("atmos.jld2", "ζ")
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
=#
