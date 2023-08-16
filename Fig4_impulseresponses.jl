using CairoMakie, PointProcessFilters, Meshes, ProgressMeter
import CairoMakie: Point2f, Circle, Rect
import PointProcessFilters: CenteredBox, CenteredBall
import Meshes: Point, Box

## set parameters
r1 = 5
r2 = 10
rsmall = (r2 - r1) / 2
defaults = (color = (:gray), strokewidth = 0)

## make figure
figure = Figure(resolution = 300 .* (1, 2 / 4), figure_padding = 1)
ax = [Axis(figure[1, i], aspect = 1, limits = (r2 + 1) .* (-1, 1, -1, 1)) for i in 1:4]
hidedecorations!.(ax)
hidespines!.(ax)
poly!(ax[1], Circle(Point2f(0, 0), r2); defaults...)
poly!(ax[1], Circle(Point2f(0, 0), r1); color = :white)

poly!(ax[2], Circle(Point2f(5, r1 + rsmall), rsmall); defaults...)
poly!(ax[2], Circle(Point2f(-5, -r1 - rsmall), rsmall); defaults...)

poly!(ax[3], Rect(-r2, -r2, 2r2, 2r2); defaults...)
poly!(ax[3], Rect(-r1, -r1, 2r1, 2r1); color = :white)

poly!(ax[4], Rect(5 - rsmall, r1, 2rsmall, 2rsmall); defaults...)
poly!(ax[4], Rect(-5 + rsmall, -r1, -2rsmall, -2rsmall); defaults...)

regions = (CenteredBall{2}(r2) - CenteredBall{2}(r1),
           Ball(Point(5, r1 + rsmall), rsmall),
           CenteredBox(Point(r2, r2)) - CenteredBox(Point(r1, r1)),
           Box(Point(5 - rsmall, r1), Point(5 + rsmall, r1 + 2rsmall)))
h = impulse_response.(RegionPassFilter.(regions))
x = range(-0.5, 0.5, 100)
h_evaluated = @showprogress "evaluating ir... " [Makie.pseudolog10.(h[i].(x, x'))
                                                 for i in 1:4]

ax = [Axis(figure[2, i], aspect = 1) for i in 1:4]
hidedecorations!.(ax)
hidespines!.(ax)
for i in 1:4
    heatmap!(ax[i], x, x, h_evaluated[i], colormap = :tempo,
             colorrange = extrema(stack(h_evaluated)), interpolate = true)
end

save("fig/fig4/transfer_functions.pdf", figure)
figure
