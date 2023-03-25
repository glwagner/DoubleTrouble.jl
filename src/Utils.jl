module Utils

using Oceananigans.Grids: minimum_xspacing, minimum_yspacing
using Oceananigans.Fields: Center

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

end # module Utils
