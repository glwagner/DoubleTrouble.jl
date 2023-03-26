module Utils

using Oceananigans
using Oceananigans.Grids: minimum_xspacing, minimum_yspacing
using Oceananigans.Operators: ζ₃ᶠᶠᶜ

const c = Center()

function split_explicit_substeps(Δt, grid; gravitational_acceleration=9.81, cfl = 0.7)
    wave_speed = sqrt(gravitational_acceleration * grid.Lz)
    min_Δx = minimum_xspacing(grid, c, c, c)
    min_Δy = minimum_yspacing(grid, c, c, c)
    Δ = 1 / sqrt(1 / min_Δx^2 + 1 / min_Δy^2)
    minimum_substeps = 10
    stable_substeps = ceil(Int, 2 * Δt / (cfl / wave_speed * Δ))
    return max(stable_substeps, minimum_substeps)
end

const OceananigansModel = Union{NonhydrostaticModel, HydrostaticFreeSurfaceModel}

function relative_vertical_vorticity(model::OceananigansModel)
    u, v, w = model.velocities
    grid = model.grid
    return KernelFunctionOperation{Face, Face, Center}(ζ₃ᶠᶠᶜ, grid, u, v)
end

function plane_simulation(arch = CPU();
                          filename,
                          size,
                          x, y,
                          z = (0, 1),
                          interface = :bottom,
                          tracers = (),
                          advection = WENO(),
                          Δt = 1)

    length(size) == 2 && (size = (size..., 1))
    grid = RectilinearGrid(arch; size, x, y, z, topology=(Periodic, Periodic, Bounded))

    # Momentum fluxes    
    boundary_conditions = Dict()
    k = interface == :bottom ? 1 : grid.Nz
    i = tuple(interface)

    τx = XFaceField(grid)
    τy = YFaceField(grid)
    τxi = interior(τx, :, :, k)
    τyi = interior(τy, :, :, k)

    ukw = NamedTuple{i}(tuple(FluxBoundaryCondition(τxi)))
    vkw = NamedTuple{i}(tuple(FluxBoundaryCondition(τyi)))
    boundary_conditions[:u] = FieldBoundaryConditions(; ukw...)
    boundary_conditions[:v] = FieldBoundaryConditions(; vkw...)

    # Tracer fluxes
    for name in tracers
        qc = CenterField(grid)
        qci = interior(qc, :, :, k)
        ckw = NamedTuple{i}(tuple(FluxBoundaryCondition(qci)))
        boundary_conditions[name] = FieldBoundaryConditions(; ckw...)
    end

    boundary_conditions = NamedTuple(name => boundary_conditions[name]
                                     for name in keys(boundary_conditions))

    model = NonhydrostaticModel(; grid, advection, tracers, boundary_conditions)
    simulation = Simulation(model; Δt=Δt, stop_time=1)
    simulation.stop_time = Inf
    wizard = TimeStepWizard(cfl=0.2, max_change=1.1)
    simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(10))

    outputs = merge(model.velocities, model.tracers, (; ζ = relative_vertical_vorticity(model)))
    simulation.output_writers[:jld2] = JLD2OutputWriter(model, outputs; filename,
                                                        schedule = IterationInterval(10),
                                                        overwrite_existing = true)

    return simulation 
end

end # module Utils
