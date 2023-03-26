module Components

struct AtmosphereOceanIce{O, A, I}
    atmosphere :: A
    ocean :: O
    ice :: I
end

struct AtmosphereOceanLand{O, A, L}
    atmosphere :: A
    ocean :: O
    land :: L
end

struct AtmosphereOcean{O, A}
    atmosphere :: A
    ocean :: O
end

#=
function set!(state::AtmosphereOcean, coupled_model::CoupledAtmosphereOceanModel)
    for name in keys(state.atmosphere)
        atmos_fields = fields(coupled_model.components.atmosphere.model)
        set!(state.atmosphere[name], atmos_fields[name])
    end

    for name in keys(state.ocean)
        ocean_fields = fields(coupled_model.components.ocean.model)
        set!(state.atmosphere[name], ocean_fields[name])
    end

    return state
end
=#

end # module Components

