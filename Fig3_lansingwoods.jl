using PointProcessFilters, CairoMakie, Meshes, ProgressMeter
import Meshes: Point, PointSet, Box
using RCall
include("src/iso_pgram.jl")
import .IsoPP: periodogram_iso

@info "Begin analysis of spatstat data"

R"""
    library(spatstat)
    treedata = lansing
"""
@rget treedata

lansing_points = Dict([s=>PointSet(unique(Point.(treedata[:x],treedata[:y])[treedata[:marks].==s])) for s in unique(treedata[:marks])])
region = Box(Point(0,0), Point(1,1))

species = ["hickory", "maple"]
λ = Dict([s=>length(lansing_points[s]) / measure(region) for s in species])
ppfilt = Dict([s=>centeredcirclepass(lansing_points[s], 6) for s in species])

lags = range(0,1,100)
pp_w = Dict(@showprogress [s=>ppfilt[s].(lags,lags').-λ[s] for s in species])

##
pgram = Dict([s=>periodogram_iso(lansing_points[s], region, 1:20) for s in species])

folder = "fig/fig3/"
clims = extrema(stack([pp_w[sp] for sp in species]))
for sp in species
    figure_p = Figure(resolution = (100,100), fontsize=8, figure_padding=1)
    ax_p = Axis(figure_p[1,1], limits=(0,1,0,1))
    hidedecorations!(ax_p)
    viz!(ax_p, lansing_points[sp], pointsize = 2, color = :black)
    save(folder * "pts_$sp.pdf", figure_p)

    figure_p = Figure(resolution = (100,60), fontsize=8, figure_padding=1)
    ax_p = Axis(figure_p[1,1], xticks = [0], yticks=[0], xgridcolor=:black, ygridcolor=:black)
    ylims!(ax_p, (-1,15))
    xlims!(ax_p, (-1,21))
    band!(ax_p, [0,6], fill(-1,2), fill(15,2), color = (:lightgrey,0.5))
    vlines!(ax_p, [0,6], color = :black, linestyle=:dash, linewidth=1/2)
    lines!(ax_p, pgram[sp].t, pgram[sp].sdf./pgram[sp].lambda.-1, color = :black)
    hidedecorations!(ax_p, grid=false)
    hidespines!(ax_p)
    save(folder * "sdf_$sp.pdf", figure_p)

    figure_p = Figure(resolution = (100,100), fontsize=8, figure_padding=1)
    ax_p = Axis(figure_p[1,1])
    hidedecorations!(ax_p)
    cf=contourf!(ax_p, lags, lags, pp_w[sp], colorrange = clims, colormap=:tempo)
    save(folder * "filt_$sp.pdf", figure_p)
end

@info "finished plotting."