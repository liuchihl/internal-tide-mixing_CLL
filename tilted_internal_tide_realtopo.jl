using Printf
using Oceananigans
using Oceananigans.Units
using Oceananigans.ImmersedBoundaries: ImmersedBoundaryGrid, GridFittedBoundary
using LinearAlgebra
using MAT
using Statistics
using Oceanostics
using Oceanostics.TKEBudgetTerms: BuoyancyProductionTerm

suffix = "0.5days"

# horizontal average 
# function terrain_following_horizontal_average(data::AbstractArray{<:Real})
#     return mean(data, dims=(1, 2))
# end


Nx = 10 #250 500 1000
Ny = 1 #500 1000 2000
Nz = 10

tᶠ = 0.5days # simulation run time
Δtᵒ = 60minutes # interval for saving output
H = 4.926kilometers # 6.e3 # vertical extent
Lx = 15kilometers
Ly = 30kilometers


## Create grid
# Creates a vertical grid with near-constant spacing `refinement * Lz / Nz` near the bottom:
# "Warped" coordinate
kwarp(k, N) = (N + 1 - k) / N
# Linear near-surface generator
ζ(k, N, refinement) = 1 + (kwarp(k, N) - 1) / refinement
# Bottom-intensified stretching function
Σ(k, N, stretching) = (1 - exp(-stretching * kwarp(k, N))) / (1 - exp(-stretching))

# Generating function
z_faces(k) = - H * (ζ(k, Nz, 1.8) * Σ(k, Nz, 10) - 1)

# using GLMakie

# lines(zspacings(grid, Center()), znodes(grid, Center()),
#       axis = (ylabel = "Depth (m)",
#               xlabel = "Vertical spacing (m)"))
# scatter!(zspacings(grid, Center()), znodes(grid, Center()), color=:red, markersize=5)


grid = RectilinearGrid(size=(Nx, Ny, Nz), 
        x = (0, Lx),
        y = (0, Ly), 
        z = z_faces,
        halo = (4,4,4),
        topology = (Periodic, Periodic, Bounded)
)
# yᶜ = ynodes(grid, Center())
# Δyᶜ = yspacings(grid, Center())

# load topography 
file = matopen("/Users/chihlunliu/Library/CloudStorage/OneDrive-SharedLibraries-UCIrvine/Chih-Lun - Documents/UC Irvine/research/topo.mat")
z_topo = read(file, "z_noslope_periodic") 
x_topo = read(file, "x_domain")
y_topo = read(file, "y_domain")
# grids has to be evenly spaced
x_topo_lin = range(x_topo[1],x_topo[end],size(z_topo,1))
y_topo_lin = range(y_topo[1],y_topo[end],size(z_topo,2))
close(file)
# high-resolution grids
x_interp = range(x_topo[1],x_topo[end], length=Nx)
# Ny=2Nx
y_interp = range(y_topo[1],y_topo[end], length=2Nx)

using Interpolations

# Interpolation object (caches coefficients and such)
itp = LinearInterpolation((x_topo_lin, y_topo_lin), z_topo)
# Interpolate
z_interp = [itp(x_topo_lin, y_topo_lin) for x_topo_lin in x_interp, y_topo_lin in y_interp]
z_interp = z_interp.-minimum(z_interp)
# heatmap(x_interp, y_interp, z_interp'; color = :balance, xlabel = "x", ylabel = "z", aspect_ratio = :equal, xlim=(0,Lx))
# z_ex = ones(Nx,Ny).*-2000
# Create immersed boundary grid

grid_real = ImmersedBoundaryGrid(grid, GridFittedBottom(z_interp[:,Nx÷2]))
velocity_bcs = FieldBoundaryConditions(immersed=ValueBoundaryCondition(0.0));

## creating terrain-aligned horizontal average
# center z grid
zc = znodes(grid, Center())
# find the grid that is above z_interp at x-y plane
inx = zeros(Nx,Ny)  # Preallocate inx array to store the indices
# create an array of indices that captures the frist element above the topography
# for i in 1:Nx
#     for j in 1:Ny
# inx[i,j] = findfirst(x -> x > z_interp[i,j], zc)
#     end
# end



# Environmental parameters
N = 1.e-3 # Brunt-Väisälä buoyancy frequency
f₀ = 0.53e-4 # Coriolis frequency
θ = 2.e-3 # tilting of domain in (x,z) plane, in radians [for small slopes tan(θ)~θ]
ĝ = [sin(θ), 0, cos(θ)] # vertical (gravity-oriented) unit vector in rotated coordinates

# Tidal forcing
U₀ = 0.025
ω₀ = 1.4e-4
u_tidal_forcing(x, y, z, t) = U₀*ω₀*sin(ω₀*t)

# IC such that flow is in phase with predicted linear response, but otherwise quiescent
Uᵣ = U₀ * ω₀^2/(ω₀^2 - f₀^2 - (N*sin(θ))^2) # quasi-resonant linear barotropic response
uᵢ(x, y, z) = -Uᵣ
vᵢ(x, y, z) = 0.
bᵢ(x, y, z) = 1e-9*rand() # seed infinitesimal perturbations in buoyancy field

s = sqrt((ω₀^2-f₀^2)/(N^2-ω₀^2))
# γ = h*π/(s*6kilometers)
# print("Steepness parameter of γ=",round(γ, digits=3))
-
# Rotate gravity vector
buoyancy = Buoyancy(model = BuoyancyTracer(), gravity_unit_vector = -ĝ)
coriolis = ConstantCartesianCoriolis(f = f₀, rotation_axis = ĝ)

# Linear background stratification (in ẑ)
@inline ẑ(x, z, ĝ) = x*ĝ[1] .+ z*ĝ[3]
@inline constant_stratification(x, y, z, t, p) = p.N² * ẑ(x, z, p.ĝ)
B̄_field = BackgroundField(constant_stratification, parameters=(; ĝ, N² = N^2))

model = NonhydrostaticModel(
    grid = grid_real,
    advection = WENO(),
    buoyancy = buoyancy,
    coriolis = coriolis,
    boundary_conditions=(u=velocity_bcs, v=velocity_bcs, w=velocity_bcs),
    forcing = (u = u_tidal_forcing,),
    closure = ScalarDiffusivity(; ν=1e-4, κ=1e-4),
    tracers = :b,
    timestepper = :RungeKutta3,
    background_fields = (; b=B̄_field),
)

set!(model, b=bᵢ, u=uᵢ, v=vᵢ)

## Configure simulation
# Δt = (1/N)*0.03
# We begin by setting the initial time step conservatively, based on the smallest grid size of our domain.
Δt = 0.5 * minimum_zspacing(grid) / Uᵣ
simulation = Simulation(model, Δt = Δt, stop_time = tᶠ)

# # The `TimeStepWizard` manages the time-step adaptively, keeping the Courant-Freidrichs-Lewy
# # (CFL) number close to `0.5` while ensuring the time-step does not increase beyond the
# # maximum allowable value for numerical stability.

wizard = TimeStepWizard(cfl=0.5, diffusive_cfl=0.2)
simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(10))


## Diagnostics
b = model.tracers.b
B̄ = model.background_fields.tracers.b
B = B̄ + b # total buoyancy field

u, v, w = model.velocities
û = @at (Face, Center, Center) u*ĝ[3] - w*ĝ[1] # true zonal velocity
ŵ = @at (Center, Center, Face) w*ĝ[3] + u*ĝ[1] # true vertical velocity

ν = model.closure.ν
κ = model.closure.κ
ε = Field(ν*(∂x(u)^2 + ∂x(v)^2 + ∂x(w)^2 + ∂y(u)^2 + ∂y(v)^2 + ∂y(w)^2 + ∂z(u)^2 + ∂z(v)^2 + ∂z(w)^2))
χ = @at (Center, Center, Center) κ[1] * (∂x(b)^2 + ∂z(b)^2)

KE = KineticEnergy(model)
ε_oceanostics = KineticEnergyDissipationRate(model)
wb = BuoyancyProductionTerm(model)
∫KE = Integral(KE)
∫ε_oceanostics = Integral(ε_oceanostics)
∫wb = Integral(wb)
# ke = @at (Center, Center, Center) 1/2 * (û^2 + ŵ^2)
# pe = @at (Center,Center,Center) -b * model.grid.zᵃᵃᶜ

# Terrain-following horizontal averages
# ε_avg = zeros(Nz-Int64(maximum(inx))+1); # includes the surface grid to the grid right above the topography
# χ_avg = zeros(Nz-Int64(maximum(inx))+1);
# ke_avg = zeros(Nz-Int64(maximum(inx))+1);
# pe_avg = zeros(Nz-Int64(maximum(inx))+1);

# for k in 1:Nz-Int64(maximum(inx))+1
# temp_ε = 0.0; temp_χ = 0.0;  temp_ke = 0.0; temp_pe = 0.0;
# for i in 1:Nx
#     for j in 1:Ny
# temp_ε = temp_ε+ε[i,j,Int64(inx[i,j])+k-1]
# temp_χ = temp_χ+χ[i,j,Int64(inx[i,j])+k-1]
# temp_ke = temp_ke+ke[i,j,Int64(inx[i,j])+k-1]
# temp_pe = temp_pe+pe[i,j,Int64(inx[i,j])+k-1]
#     end
# end
# ε_avg[k] = temp_ε / (Nx*Ny)
# χ_avg[k] = temp_χ / (Nx*Ny)
# ke_avg[k] = temp_ke / (Nx*Ny)
# pe_avg[k] = temp_pe / (Nx*Ny)
# end
# ε_avg_func(m) = m.ε_avg  # Define ε_avg function
# χ_avg_func(m) = m.χ_avg  # Define χ_avg function
# ε_avg_field = AbstractVector{ε_avg_func(m)}
# χ_avg_field = AbstractVector{χ_avg_func}



# output_data = NamedTuple(:ε_avg => ε_avg_func, :χ_avg => χ_avg_func)

# ε_avg_func(m) = m.ε_avg  # Define ε_avg function
# χ_avg_func(m) = m.χ_avg  # Define χ_avg function
# output_data = NamedTuple(:ε_avg => (m) -> ε_avg(m), :χ_avg => (m) -> χ_avg(m))

custom_diags = (B=B, uhat=û, what=ŵ, χ=χ, ε=ε, KE, ε_oceanostics, wb, ∫KE, ∫ε_oceanostics, ∫wb)
all_diags = merge(model.velocities, model.tracers, custom_diags)
# avg_diags = (ε_avg_field=ε_avg_field, χ_avg_field=χ_avg_field)

# ε_avg_output(model,location=(center,center,center)) = ε_avg
# χ_avg_output(model,location=(center,center,center)) = χ_avg
# ke_avg_output(model,location=(center,center,center)) = ke_avg
# pe_avg_output(model,location=(center,center,center)) = pe_avg
# avg_diags = NamedTuple(:ε_avg => ε_avg_output(model), :χ_avg => χ_avg_output(model), :ke_avg => ke_avg_output(model), :pe_avg => pe_avg_output(model))
# avg_diags = NamedTuple(:ε_avg => ε_avg_output, :χ_avg => χ_avg_output, :ke_avg => ke_avg_output, :pe_avg => pe_avg_output)

fname = string("internal_tide_", suffix,"-theta=",string(θ),"2Drealtopo_coarse")

# JLD2OutputWriter  
simulation.output_writers[:checkpointer] = Checkpointer(
                                        model,
                                        schedule=TimeInterval(Δtᵒ*2),
                                        dir="output",
                                        prefix=string(fname, "_checkpoint"),
                                        cleanup=true)

simulation.output_writers[:fields] = JLD2OutputWriter(model, all_diags,
                                        schedule = TimeInterval(Δtᵒ*2),
                                        filename = string("output/", fname, "_state.jld2"),
                                        overwrite_existing = true)

simulation.output_writers[:section_snapshots] = JLD2OutputWriter(model, (; b, χ),
                                        schedule = TimeInterval(Δtᵒ),
                                        indices = (:,1,:),
                                        filename = string("output/", fname, "_section_snapshots.jld2"),
                                        overwrite_existing = true)

# simulation.output_writers[:zonal_time_means] = JLD2OutputWriter(model, (; ε),
#                                         schedule = AveragedTimeInterval(Δtᵒ÷2, window=Δtᵒ÷2),
#                                         filename = string("output/", fname, "_zonal_time_means.jld2"),
#                                         overwrite_existing = true)

# simulation.output_writers[:TF_horizontal_average] = JLD2OutputWriter(model, avg_diags;
#                                         schedule = AveragedTimeInterval(Δtᵒ÷2, window=Δtᵒ÷2),
#                                         filename = string("output/", fname, "_TF_horizontal_average.jld2"),
#                                         overwrite_existing = true)

progress_message(s) = @info @sprintf("[%.2f%%], iteration: %d, time: %.3f, max|w|: %.2e, 
                            advective CFL: %.2e, diffusive CFL: %.2e\n",
                            100 * s.model.clock.time / s.stop_time, s.model.clock.iteration,
                            s.model.clock.time, maximum(abs, model.velocities.w),
                            AdvectiveCFL(s.Δt)(s.model), DiffusiveCFL(s.Δt)(s.model))
simulation.callbacks[:progress] = Callback(progress_message, TimeInterval(Δtᵒ))

## Running the simulation!
# run!(simulation)
## Running the simulation!
run!(simulation)

@info """
    Simulation complete.
    Output: $(abspath(simulation.output_writers[:fields].filepath))
"""

