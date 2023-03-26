module DoubleTrouble

export CoupledModel

using Oceananigans.Grids: AbstractGrid, Center, Face
using Oceananigans.Utils: prettysummary
using Oceananigans.Fields: CenterField, XFaceField, YFaceField, interior, set!
using Oceananigans.Models: AbstractModel, NonhydrostaticModel, HydrostaticFreeSurfaceModel
using Oceananigans.TimeSteppers: Clock, tick!
using Oceananigans.AbstractOperations: KernelFunctionOperation
using Oceananigans.Operators: ζ₃ᶠᶠᶜ
using Oceananigans.Simulations: run!

# Simulations interface
import Oceananigans: prognostic_fields
import Oceananigans.TimeSteppers: time_step!, update_state!
import Oceananigans.Models: initialize_model!
import Oceananigans.Simulations: reset!

include("Utils.jl")
include("Fluxes.jl")
include("Components.jl")

using .Components: AtmosphereOcean

reset!(::Nothing) = nothing # silly Oceananigans # ⟨⟨ eyes ⟩⟩

#####
##### CoupledModel
#####

mutable struct CoupledModel{Comp, Flux, Surf, FT}
    grid :: Nothing # just a place holder
    timestepper :: Nothing # just a place holder
    components :: Comp
    surfaces :: Surf
    fluxes :: Flux
    clock :: Clock{FT}
end

const CoupledAtmosphereOceanModel = CoupledModel{<:AtmosphereOcean}

function Base.show(io::IO, m::CoupledAtmosphereOceanModel)

    Nfluxes = length(m.fluxes)
    if Nfluxes > 0
        item  = "├── "
        space = "│    "
    else
        item  = "└── "
        space = "     "
    end

    FT = typeof(m.clock.time)
    print(io, "CoupledAtmosphereOceanModel{$FT}:", '\n')
    print(io, "├── atmosphere: ", prettysummary(m.components.atmosphere.model), '\n')
    print(io, "│   └── grid: ", prettysummary(m.components.atmosphere.model.grid), '\n')
    print(io, item, "ocean: ", prettysummary(m.components.ocean.model), '\n')
    print(io, space, "└── grid: ", prettysummary(m.components.ocean.model.grid))

    if Nfluxes > 0
        print(io, '\n')
        if Nfluxes > 1
            print(io, "└── ", Nfluxes, " fluxes: ", '\n')
        else
            print(io, "└── ", Nfluxes, " flux: ", '\n')
        end
    end

    for (n, flux) in enumerate(m.fluxes)
        if n < Nfluxes
            print(io, "    ├── ", prettysummary(flux), '\n')
        else
            print(io, "    └── ", prettysummary(flux))
        end
    end
end

tupleit(t::Tuple) = t

tupleit(t) =
    try
        Tuple(t)
    catch
        tuple(t)
    end

function CoupledAtmosphereOceanModel(; ocean, atmosphere, fluxes)
    components = AtmosphereOcean(atmosphere, ocean)
    fluxes = tupleit(fluxes)

    #=
    # Placeholder surface state?
    oNx, oNy, oNz = size(ocean.model.grid)
    ocean_u = view(ocean.model.velocities.u, :, :, oNz)
    ocean_v = view(ocean.model.velocities.v, :, :, oNz)
    ocean_surface = (u=ocean_u, v=ocean_v)

    atmos_u = view(atmosphere.model.velocities.u, :, :, 1)
    atmos_v = view(atmosphere.model.velocities.v, :, :, 1)
    atmos_surface = (u=atmos_u, v=atmos_v)

    surfaces = AtmosphereOcean(atmos_surface, ocean_surface)
    =#
    surfaces = nothing

    clock = Clock(time=zero(ocean.model.grid))

    return CoupledModel(nothing, nothing, components, surfaces, fluxes, clock)
end

# Oceananigans.Simulations interface
prognostic_fields(m::CoupledModel) = prognostic_fields(m.components.atmos.model)
initialize_model!(::CoupledModel) = nothing

function update_state!(coupled_model::CoupledModel)

    # Regridding, computation of surface state, etc happens here.
    # Not always needed.

    for flux in coupled_model.fluxes
        compute_flux!(flux, coupled_model)
    end

    # Coupler "push": distribute fluxes to components
    set!(coupled_model.components, coupled_model.fluxes)

    return nothing
end
                       
function time_step!(model::CoupledAtmosphereOceanModel, Δt; callbacks=nothing)
    simulations = (model.components.ocean, model.components.atmos)

    time = model.clock.time
    stop_time = model.clock.time + Δt

    # Run components asynchronously
    asyncmap(simulations) do simulation
        simulation.verbose = false
        simulation.model.clock.time = time
        simulation.stop_time = stop_time
        run!(simulation)
    end

    update_state!(model)
    tick!(model.clock, Δt)

    return nothing
end

end # module DoubleTrouble
