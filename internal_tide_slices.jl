using Printf
using Oceananigans
using Oceananigans.Units
using Oceananigans.TurbulenceClosures: VerticallyImplicitTimeDiscretization
using Oceananigans.ImmersedBoundaries: ImmersedBoundaryGrid, GridFittedBoundary
# using CairoMakie
using JLD2
using Plots

function nice_divergent_levels(c, clim; nlevels=20)
    levels = range(-clim, stop=clim, length=nlevels)
    cmax = maximum(abs, c)
    clim < cmax && (levels = vcat([-cmax], levels, [cmax]))
    return (-clim, clim), levels
end

function nan_solid(x, z, u, mask)
    Nx, Nz = size(u)
    x2 = reshape(x, (Nx, 1))
    z2 = reshape(z, (1, Nz))
    u[mask] .= NaN
    return u
end

saved_output_filename = "output/internal_tide_0.5days-theta=0.0022Drealtopo_coarse_zonal_time_means.jld2"
# saved_output_filename = "internal_tide_sinusoidal_20-days.jld2"

file = jldopen(saved_output_filename)

ig = file["serialized/grid"]
ug = ig.underlying_grid
ĝ = file["serialized/buoyancy"].gravity_unit_vector

û = FieldTimeSeries(saved_output_filename, "uhat")
ŵ = FieldTimeSeries(saved_output_filename, "what")
# ε = FieldTimeSeries(saved_output_filename, "ε")
χ = FieldTimeSeries(saved_output_filename, "χ")
B = FieldTimeSeries(saved_output_filename, "B")
t = B.times

xu, yu, zu = nodes(û[1])
xw, yw, zw = nodes(ŵ[1])
xb, yb, zb = nodes(B[1])

# for n in 1:length(t)
#     ε[n].data .= log10.(@. ifelse(ε[n].data > 0, ε[n].data, 1.e-100)) # transform ε to logspace
# end

# dBdz = Field(∂z(B))
# heatmap!(xu,zu,dBdz[:,1,:,500])
# heatmap(xb, zb, dBdz[:,1,:,500]'; color = :thermal, xlabel = "x", ylabel = "z", aspect_ratio = :equal); 


# making some plots
u_ic = FieldTimeSeries( saved_output_filename, "u", iterations = 0)

U₀ = 0.025
û[û .== 0] .= NaN
hmu = heatmap(xu , zu, û[1:end,1,1:end,25]',
 xlabel = "x", ylabel = "z (m)",colorrange = (-3U₀, 3U₀),
 clim=(-3U₀, 3U₀), colormap = :balance,lowclip=cgrad(:balance)[1], 
 highclip=cgrad(:balance)[end])

 ŵ[ŵ .== 0] .= NaN
hmw = heatmap(xw , zw, ŵ[1:end,1,1:end,25]', 
xlabel = "x", ylabel = "z (m)",colorrange = (-3U₀, 3U₀),
clim=(-3U₀, 3U₀))


# ε[ε .== 0] .= NaN
# hmeps = heatmap(xw , zw, log10(ε[1:end,1,1:end,25]'), 
# xlabel = "x", ylabel = "z (m)",colorrange = (-3U₀, 3U₀),
# clim=(-3U₀, 3U₀))
# plot(hmu, hmw, layout=(2,1))

χ[χ .== 0] .= NaN
hmchi = heatmap(xb , zb, log10.(χ[1:end,1,1:end,25]'), 
xlabel = "x (m)", ylabel = "z (m)",colorrange = (-20,-0.5),
clim=(-8, -2))

plot(hmu, hmw, hmchi, layout=(3,1), size = (1200, 800))
