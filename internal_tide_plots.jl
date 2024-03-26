using Printf
using Oceananigans
using Oceananigans.Units
using Oceananigans.TurbulenceClosures: VerticallyImplicitTimeDiscretization
using Oceananigans.ImmersedBoundaries: ImmersedBoundaryGrid, GridFittedBoundary
using CairoMakie
using JLD2
# using Plots

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

saved_output_filename = "output/internal_tide_3days-theta=0.0036_realtopo3D_Nx150_section_snapshots.jld2"
# saved_output_filename = "internal_tide_sinusoidal_20-days.jld2"

file = jldopen(saved_output_filename)

ig = file["serialized/grid"]
ug = ig.underlying_grid
ĝ = file["serialized/buoyancy"].gravity_unit_vector

# û = FieldTimeSeries(saved_output_filename, "uhat")
# ŵ = FieldTimeSeries(saved_output_filename, "what")
# ε = FieldTimeSeries(saved_output_filename, "ε")
# B = FieldTimeSeries(saved_output_filename, "B")
# t = B.times

wb = FieldTimeSeries(saved_output_filename, "wb")
ε_oceanostics = FieldTimeSeries(saved_output_filename, "ε_oceanostics")
KE = FieldTimeSeries(saved_output_filename, "KE")
χ = FieldTimeSeries(saved_output_filename, "χ")
zc = ug.zᵃᵃᶜ[1:end]
xc = ug.xᶜᵃᵃ[1:end]
do i [0,3:2:]
wb[i] = file["timeseries/wb/$i"]
end

xu, yu, zu = nodes(û[1])
xw, yw, zw = nodes(ŵ[1])
xb, yb, zb = nodes(B[1])

for n in 1:length(t)
    ε[n].data .= log10.(@. ifelse(ε[n].data > 0, ε[n].data, 1.e-100)) # transform ε to logspace
end

n = Observable(1)

# mask immersed boundaries
# mask = interior(ig.immersed_boundary.mask)[:,1,:]
# ûₙ = @lift nan_solid(xu, zu, interior(û[$n], :, 1, :), mask)
# εₙ = @lift nan_solid(xb, zb, interior(ε[$n], :, 1, :), mask)
# Bₙ = @lift nan_solid(xb, zb, interior(B[$n], :, 1, :), mask)
ûₙ = @lift interior(û[$n], :, 1, :)
εₙ = @lift interior(ε[$n], :, 1, :)
Bₙ = @lift interior(B[$n], :, 1, :)
begin
    fig = Figure(resolution = (1000, 1000), figure_padding=(10, 40, 10, 10))
    axis_kwargs = (xlabel = "zonal distance (x)",
                      ylabel = "elevation (z)",
                      limits = ((0, ug.Lx), (0, ug.Lz)),
                      )
    ax_u = Axis(fig[2, 1]; title = "zonal velocity (u) and equally-spaced buoyancy contours (B)", axis_kwargs...)
    ax_ε = Axis(fig[3, 1]; title = "TKE dissipation rate (ε) and equally-spaced buoyancy contours (B)", axis_kwargs...)
    # ax_dbdz = Axis(fig[3, 1]; title = "dB/dz", axis_kwargs...)
end

U₀ = 0.025
hm_u = heatmap!(ax_u, xu, zu, ûₙ,
    colorrange = (-3U₀, 3U₀), colormap = :balance,
    lowclip=cgrad(:balance)[1], highclip=cgrad(:balance)[end])
ct_u = contour!(ax_u, xb, zb, Bₙ,
    levels=0.:0.5e-4:4.e-3, linewidth=0.25, color=:black, alpha=0.2)

hm_ε = heatmap!(ax_ε, xb, zb, εₙ,
    colorrange = (-11, -8), colormap = :matter,
    lowclip=cgrad(:matter)[1], highclip=cgrad(:matter)[end])
ct_ε = contour!(ax_ε, xb, zb, Bₙ,
    levels=0.:0.5e-4:4.e-3, linewidth=0.25, color=:black, alpha=0.2)

frames = (1:length(t))
record(fig, string(split(saved_output_filename, "1.")[1], ".mp4"), frames, framerate=16) do i
    @info "Plotting frame $i of $(frames[end])..."
    n[] = i
end

# dBdz = Field(∂z(B))
# heatmap!(xu,zu,dBdz[:,1,:,500])
# heatmap(xb, zb, dBdz[:,1,:,500]'; color = :thermal, xlabel = "x", ylabel = "z", aspect_ratio = :equal); 


# # making some plots
# u_ic = FieldTimeSeries( saved_output_filename, "u", iterations = 0)

# hmu = heatmap(xu , zu, û[1:end,1,1:end,25], xlabel = "x", ylabel = "z (m)")
# hmw = heatmap!(xw , zw, ŵ[1:end,1,1:end,25], xlabel = "x", ylabel = "z (m)")

# plot(hmu, hmw, layout=(2,1))