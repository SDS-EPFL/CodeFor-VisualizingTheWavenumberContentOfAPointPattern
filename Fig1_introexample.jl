using PointProcessFilters, CairoMakie, Meshes, ProgressMeter, RCall, LinearAlgebra

## Simulation
@info "Perform the simulation..."
side_length = 200
side_padding = 50
box_x = [-side_padding, side_length + side_padding]
box_y = box_x

λ = 0.01
kappa = 0.8λ
scale = 3
mu = λ / kappa

@rput kappa scale mu box_x box_y
R"""
    library(spatstat)
    set.seed(1234)
    simulated_process = rThomas(kappa, scale, mu, win=owin(box_x, box_y))
"""
@rget simulated_process
simulated_points = PointSet(Meshes.Point.(simulated_process[:x], simulated_process[:y]))
@info "points simulated."

## Filtering
@info "Apply the filters..."
interval = [0, 0.05, 0.1, 0.15]
filtered_process = [
    centeredcirclepass(simulated_points, interval[2]),
    ringpass(simulated_points, interval[2], interval[3]),
    ringpass(simulated_points, interval[3], interval[4]),
]
filter_resolution = 100
filter_x = range(0, side_length, filter_resolution)
filter_y = filter_x
Y = @showprogress "evaluating filters... " [y.(filter_x, filter_y')
                                            for y in filtered_process]
Ylims = extrema(stack(Y))
@info "filtered processes computed."

## sdf of thomas process
function thomassdf(k; kappa, scale, mu)
    return kappa * mu * (1 + mu * exp(-scale^2 * (2π)^2 * norm(k)^2))
end

## plotting
@info "Start plotting..."
fig_res = (80, 80)

## plot the spectra
@info "- plotting spectra"
figure = Figure(resolution = fig_res .* (1, 2 / 3), figure_padding = 1)
ax = Axis(figure[1, 1], yscale = identity, xticks = [0], yticks = [0], xgridcolor = :black,
          ygridcolor = :black)
band!(ax, interval, zeros(length(interval)), fill(0.025, length(interval)),
      color = (:lightgrey, 0.5))
vlines!(interval, color = :black, linestyle = :dash, linewidth = 1 / 2)
lines!(ax, range(0, 0.2, 1000), x -> thomassdf(x, kappa = kappa, scale = scale, mu = mu),
       color = :black)

ylims!(ax, (-0.001, 0.025))
hidedecorations!(ax, grid = false)
hidespines!(ax)
resize_to_layout!(figure)

save("fig/fig1/sdf.pdf", figure)
figure

## plot the points
@info "- plotting points"
figure = Figure(resolution = fig_res, figure_padding = 1)
ax = Axis(figure[1, 1], limits = (0, side_length, 0, side_length))
viz!(ax, simulated_points, color = :black, pointsize = 2)
hidedecorations!(ax)
resize_to_layout!(figure)
save("fig/fig1/points.pdf", figure)
figure

## plot filters
@info "- plotting filters"
for i in eachindex(Y)
    figure_filt = Figure(resolution = fig_res, figure_padding = 1)
    ax_filt = Axis(figure_filt[1, 1])
    contourf!(ax_filt, filter_x, filter_y, Y[i], colorrange = Ylims, colormap = :tempo)
    hidedecorations!(ax_filt)
    save("fig/fig1/filtered_$i.pdf", figure_filt)
end

@info "Finished plotting."
