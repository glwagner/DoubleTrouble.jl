module Fluxes

using Oceananigans
using Oceananigans.Fields: AbstractField
using Oceananigans.Grids: AbstractGrid

import Oceananigans.Fields: compute!
import Oceananigans.Fields: set!

struct Flux{F, D}
    formula :: F
    field :: D
end

#####
##### Formula
#####

abstract type AbstractBulkFormula end

#####
##### Momentum flux
#####

Base.@kwdef struct AirSeaQuadraticDrag <: AbstractBulkFormula
    drag_coefficient :: Float64 = 2e-3
end

const AirSeaMomentumFlux = Flux{<:AirSeaQuadraticDrag}

Base.summary(::AirSeaMomentumFlux) = "AirSeaMomentumFlux"

struct SurfaceMomentumField{U, V}
    u :: U
    v :: V
end

function relative_air_sea_velocities(coupled_model)
    uo = coupled_model.components.ocean.model.velocities.u
    vo = coupled_model.components.ocean.model.velocities.v
    ua = coupled_model.components.atmosphere.model.velocities.u
    va = coupled_model.components.atmosphere.model.velocities.v
    Δu = uo - ua
    Δv = vo - va
    return (u=Δu, v=Δv)
end

air_density(atmos) = 1.0
seawater_density(ocean) = 1024.0

function compute_flux!(τ::AirSeaMomentumFlux, coupled_model)
    Δu, Δv = relative_air_sea_surface_velocities(coupled_model)
    ΔU = sqrt(Δu^2 + Δv^2)

    Cᴰ = τ.formula.drag_coefficient
    ρ_ocean = seawater_density(coupled_model.ocean)
    ρ_atmos = air_density(coupled_model.atmos)

    τ.field.u .= Cᴰ * ρ_atmos * Δu * ΔU
    τ.field.v .= Cᴰ * ρ_atmos * Δv * ΔU

    return nothing
end

function set!(components::AtmosphereOcean, τ::AirSeaMomentumFlux)

    # TODO: these should be separate functions
    # So model components will "support" `set!`
    #
    # Questions:
    #   - How to identify components? "Ocean", "Atmosphere", etc.

    ocean = components.ocean
    ρ_ocean = seawater_density(coupled_model.ocean)
    uo, vo, wo = ocean.model.velocities
    τox = uo.boundary_conditions.top.condition
    τoy = vo.boundary_conditions.top.condition
    interior(τox) .= interior(τ.field.u) ./ ρ_ocean
    interior(τoy) .= interior(τ.field.v) ./ ρ_ocean

    atmos = components.atmosphere
    ρ_atmos = air_density(coupled_model.atmos)
    ua, va, wa = atmos.model.velocities
    τax = ua.boundary_conditions.bottom.condition
    τay = va.boundary_conditions.bottom.condition
    interior(τax) .= interior(τ.field.u) ./ ρ_atmos
    interior(τay) .= interior(τ.field.v) ./ ρ_atmos

    return nothing
end

function Flux(drag_formula::AirSeaQuadraticDrag, grid::AbstractGrid)
    τx = XFaceField(grid)
    τy = YFaceField(grid)
    data = SurfaceMomentumField(τx, τy)
    return Flux(drag_formula, data)
end

#####
##### Tracer flux
#####

struct AirSeaTracerTransfer <: AbstractBulkFormula
    tracer_name :: Symbol
    transfer_coefficient :: Float64
end

const AirSeaTracerFlux = Flux{<:AirSeaTracerTransfer}
Base.summary(q::AirSeaTracerFlux) = string("AirSeaTracerFlux of ", q.formula.tracer_name)

function Flux(transfer_formula::AirSeaTracerTransfer, grid::AbstractGrid)
    data = CenterField(grid)
    return Flux(transfer_formula, data)
end

function compute_flux!(q::AirSeaTracerFlux, coupled_model)
    name = q.formula.tracer_name
    co = coupled_model.components.ocean.tracers[name]
    ca = coupled_model.components.atmosphere.tracers[name]
    Δc = co - ca

    ρ_ocean = seawater_density(coupled_model.ocean)
    ρ_atmos = air_density(coupled_model.atmos)

    Δu, Δv = relative_air_sea_surface_velocities(coupled_model)
    ΔU = sqrt(Δu^2 + Δv^2)

    C = q.formula.transfer_coefficient
    q.field .= C * ρ_atmos * Δc * ΔU

    return nothing
end

function set!(components::AtmosphereOcean, q::AirSeaTracerFlux)

    # TODO: these should be separate functions
    # So model components will "support" `set!`
    #
    # Questions:
    #   - How to identify components? "Ocean", "Atmosphere", etc.

    name = q.formula.tracer_name

    ocean = components.ocean
    ρ_ocean = seawater_density(coupled_model.ocean)
    co = ocean.model.tracers[name]
    qo = co.boundary_conditions.top.condition
    interior(qo) .= interior(q.field) ./ ρ_ocean

    atmos = components.atmosphere
    ρ_atmos = air_density(coupled_model.atmos)
    ca = atmos.model.tracers[name]
    qa = ca.boundary_conditions.bottom.condition
    interior(qa) .= interior(q.field) ./ ρ_atmos
    
    return nothing
end

end # module BulkFormulae

