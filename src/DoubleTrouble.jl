module DoubleTrouble

export CoupledModel

using Oceananigans.Grids: AbstractGrid, Center, Face
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

reset!(::Nothing) = nothing # silly Oceananigans # ⟨⟨ eyes ⟩⟩

const OceananigansModel = Union{NonhydrostaticModel, HydrostaticFreeSurfaceModel}

function relative_vertical_vorticity(model::OceananigansModel)
    u, v, w = model.velocities
    grid = model.grid
    return KernelFunctionOperation{Face, Face, Center}(ζ₃ᶠᶠᶜ, grid, u, v)
end

"""
    AtmosphereOcean

Generic type for holding a tuple of (atmosphere, ocean) types: simulations, states, fluxes.
"""
struct AtmosphereOcean{O, A}
    atmos :: A
    ocean :: O
end

function AtmosphereOcean(atmos_grid::AbstractGrid, ocean_grid::AbstractGrid=atmos_grid)
    atmos = (u = XFaceField(atmos_grid),
             v = YFaceField(atmos_grid))

    ocean = (u = XFaceField(ocean_grid),
             v = YFaceField(ocean_grid))

    return AtmosphereOcean(atmos, ocean)
end

#####
##### SurfaceFluxes: calculates bulk formula
#####

struct SurfaceFluxes{F, B, S, G}
    grid :: G
    states :: S
    fluxes :: F
    bulk_formulae :: B
end

const AtmosphereOceanSurfaceFluxes = SurfaceFluxes{<:AtmosphereOcean}

Base.@kwdef struct SimpleBulkFormulae
    drag_coefficient :: Float64 = 2e-3
    ocean_density :: Float64 = 1024.0
    atmos_density :: Float64 = 1.0
end

function update_state!(sf::AtmosphereOceanSurfaceFluxes, components::AtmosphereOcean)
    ocean_model = components.ocean.model
    atmos_model = components.atmos.model

    oNx, oNy, oNz = size(ocean_model.grid)
    ocean_surface_u = interior(ocean_model.velocities.u, :, :, oNz)
    ocean_surface_v = interior(atmos_model.velocities.v, :, :, oNz)

    aNx, aNy, aNz = size(atmos_model.grid)
    atmos_surface_u = interior(atmos_model.velocities.u, :, :, 1)
    atmos_surface_v = interior(atmos_model.velocities.v, :, :, 1)

    set!(sf.states.ocean.u, ocean_surface_u)
    set!(sf.states.ocean.v, ocean_surface_v)

    set!(sf.states.atmos.u, atmos_surface_u)
    set!(sf.states.atmos.v, atmos_surface_v)

    return nothing
end

function compute_surface_fluxes!(fluxes::AtmosphereOcean,
                                 bf::SimpleBulkFormulae,
                                 states::AtmosphereOcean)

    uo = interior(states.ocean.u)
    vo = interior(states.ocean.v)
    ua = interior(states.atmos.v)
    va = interior(states.atmos.v)

    # These are arrays for now
    # but should be Field in future
    τax = fluxes.atmos.u
    τay = fluxes.atmos.v
    τox = fluxes.ocean.u
    τoy = fluxes.ocean.u

    ρa = bf.atmos_density
    ρo = bf.ocean_density
    Cᴰ = bf.drag_coefficient

    @. τax = Cᴰ * (uo - ua)^2
    @. τay = Cᴰ * (vo - va)^2

    # Update the surface fluxes in component models
    @. τox = ρa / ρo * τax
    @. τoy = ρa / ρo * τay

    return nothing
end

function top_momentum_fluxes(model::OceananigansModel)
    velocities = model.velocities
    u_flux = velocities.u.boundary_conditions.top.condition     
    v_flux = velocities.v.boundary_conditions.top.condition     
    return (u=u_flux, v=v_flux)
end

function bottom_momentum_fluxes(model::OceananigansModel)
    velocities = model.velocities
    u_flux = velocities.u.boundary_conditions.bottom.condition     
    v_flux = velocities.v.boundary_conditions.bottom.condition     
    return (u=u_flux, v=v_flux)
end

"""
    SurfaceFluxes(surface_grid; kw...)

Return a utility for computing `fluxes` on `surface_grids` using
`bulk_formulae`, which are functions of the component `states`.
"""
function SurfaceFluxes(surface_grid, ocean, atmos;
                       bulk_formulae = SimpleBulkFormulae(),
                       states = AtmosphereOcean(surface_grid),
                       fluxes = nothing)

    if isnothing(fluxes) # extract automagically
        τ_ocean = top_momentum_fluxes(ocean.model)
        τ_atmos = bottom_momentum_fluxes(atmos.model)
        fluxes = AtmosphereOcean(τ_atmos, τ_ocean)
    end

    return SurfaceFluxes(surface_grid, states, fluxes, bulk_formulae)
end

#####
##### CoupledModel
#####

struct CoupledModel{Comp, Flux, FT, Grid}
    grid :: Grid
    timestepper :: Nothing # just a place holder
    surface_fluxes :: Flux
    components :: Comp
    clock :: Clock{FT}
end

const CoupledAtmosphereOceanModel = CoupledModel{<:AtmosphereOcean}

function CoupledAtmosphereOceanModel(surface_grid; ocean, atmosphere)
    surface_fluxes = SurfaceFluxes(surface_grid, ocean, atmosphere)
    components = AtmosphereOcean(atmosphere, ocean)
    clock = Clock(time=zero(surface_grid))
    return CoupledModel(surface_grid, nothing, surface_fluxes, components, clock)
end

# Oceananigans.Simulations interface
prognostic_fields(m::CoupledModel) = prognostic_fields(m.components.atmos.model)
initialize_model!(::CoupledModel) = nothing

function update_state!(coupled_model::CoupledModel)
    sf = coupled_model.surface_fluxes
    update_state!(sf, coupled_model.components)
    compute_surface_fluxes!(sf.fluxes, sf.bulk_formulae, sf.states)
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
