using CairoMakie, PointProcessFilters, ProgressMeter, GaussianRandomFields, RCall, Meshes
import Meshes: Point, PointSet, Box
include("src/iso_pgram.jl")
import .IsoPP: periodogram_iso

## set parameter
side_length = 100
side_padding = 25
box_x = [-side_padding, side_length + side_padding]
box_y = box_x
nu = 3 / 2
lengthscale = 10
sigma = sqrt(2)
lambda_noise = 0.05

## simualate the random field (R package currently broken so cant use spatstat)
import Random;
Random.seed!(1234);
cov = CovarianceFunction(2, Matern(lengthscale, nu, σ = sigma))
grid_sides = range(box_x..., 250), range(box_y..., 250)
grf = GaussianRandomField(cov, CirculantEmbedding(), grid_sides...)
intensity = exp.(GaussianRandomFields.sample(grf) .- 4)
heatmap(intensity)

##
@rput intensity box_x box_y lambda_noise
R"""
    library(spatstat)
    set.seed(1234)
    window = owin(box_x, box_y)
    intensity_image = as.im(aperm(intensity), window)
    signal_process = rpoispp(intensity_image, nsim=2)
    noise_process = rpoispp(lambda_noise, win=window, nsim=2)
"""
@rget signal_process noise_process
sim_names = collect(keys(signal_process))
signal = Dict([k => PointSet(Point.(signal_process[k][:x], signal_process[k][:y]))
               for k in sim_names])
noisy = Dict([k => PointSet([Point.(signal_process[k][:x], signal_process[k][:y]);
                             Point.(noise_process[k][:x], noise_process[k][:y])])
              for k in sim_names])

## apply filters
@info "Filtering (may take a few minutes)..."
area = diff(box_x)[1] * diff(box_y)[1]
filtered = Dict([k => centeredcirclepass(noisy[k], 0.05).(grid_sides[1], grid_sides[2]') .-
                      length(noisy[k]) / area for k in sim_names])
@info "filtering completed."

##
clims = extrema(vcat(values(filtered)...))
f = Figure()
heatmap(f[1, 1], grid_sides..., intensity, colormap = :matter)
for i in 1:2
    viz(f[2, i], signal[sim_names[i]], pointsize = 4, color = :black,
        axis = (limits = (0, 100, 0, 100),))
    viz(f[3, i], noisy[sim_names[i]], pointsize = 4, color = :black,
        axis = (limits = (0, 100, 0, 100),))
    heatmap(f[4, i], grid_sides..., filtered[sim_names[i]], colormap = :tempo,
            colorrange = clims)
end
f

region = Box(Point(box_x[1],box_y[1]), Point(box_x[2],box_y[2]))
pgram = Dict([s=>periodogram_iso(noisy[s], region, (1:20).*0.01) for s in sim_names])
λ = Dict([s=>length(noisy[s]) / area for s in sim_names])

##
@info "Making figures..."
clims = extrema(vcat(values(filtered)...))
function make_figure()
    f = Figure(resolution = (100, 100), figure_padding = 1)
    ax = Axis(f[1, 1], limits = (0, 100, 0, 100), aspect = 1)
    hidedecorations!(ax)
    return f, ax
end

figure, ax = make_figure()
contourf!(ax, grid_sides..., intensity, colormap = :matter)
save("fig/fig2/intensity.pdf", figure)

for i in 1:2
    figure_p, ax_p = make_figure()
    viz!(ax_p, signal[sim_names[i]], pointsize = 2, color = :black)
    save("fig/fig2/cp$i.pdf", figure_p)

    figure_p, ax_p = make_figure()
    viz!(ax_p, noisy[sim_names[i]], pointsize = 2, color = :black)
    save("fig/fig2/pts$i.pdf", figure_p)

    figure_p, ax_p = make_figure()
    contourf!(ax_p, grid_sides..., filtered[sim_names[i]], colormap = :tempo,
              colorrange = clims)
    save("fig/fig2/filt$i.pdf", figure_p)

    figure_p = Figure(resolution = (100,60), fontsize=8, figure_padding=1)
    ax_p = Axis(figure_p[1,1], xticks = [0], yticks=[0], xgridcolor=:black, ygridcolor=:black)
    ylims!(ax_p, (-1,21))
    band!(ax_p, [0,0.05], fill(-1,2), fill(21,2), color = (:lightgrey,0.5))
    vlines!(ax_p, [0,0.05], color = :black, linestyle=:dash, linewidth=1/2)
    lines!(ax_p, pgram[sim_names[i]].t, pgram[sim_names[i]].sdf./pgram[sim_names[i]].lambda.-1, color = :black)
    hidedecorations!(ax_p, grid=false)
    hidespines!(ax_p)
    save("fig/fig2/sdf_$i.pdf", figure_p)
end
@info "figures finished."
