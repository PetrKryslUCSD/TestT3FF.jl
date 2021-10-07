using LinearAlgebra
using FinEtools
using PGFPlotsX

objects = []

E = 200*phun("GPa");
nu = 0.3;
G = E/2/(1+nu)
a = 1.0*phun("m");
t = a/5
D = E*t^3/12/(1-nu^2)
P = 1.0
a = 1.0
r = collect(t/100000:t/10000:t)
w_b = @. P/16/pi/D * (2*r^2 * log(r/a) + (3+nu)/(1+nu)*(a^2 - r^2))
w_s = @. -P/4/pi/G/t * log(r/a) 
max_w_b = maximum(w_b)

@pgf p = PGFPlotsX.Plot(
{
color = "black",
line_width  = 0.7, 
style = "solid"
},
Coordinates([v for v in  zip(r./t, w_b./max_w_b)])
)
push!(objects, p)
push!(objects, LegendEntry("Bending only"))


@pgf p = PGFPlotsX.Plot(
{
color = "red",
line_width  = 0.7, 
style = "solid"
},
Coordinates([v for v in  zip(r./t, (w_b+w_s)./max_w_b)])
)
push!(objects, p)
push!(objects, LegendEntry("Bending plus shear"))



@pgf ax = Axis(
    {
        xlabel = "Normalized Distance [ND]",
        ylabel = "Normalized Displacement [ND]",
        # xmin = range[1],
        # xmax = range[2],
        xmode = "linear", 
        ymode = "linear",
        yminorgrids = "true",
        grid = "both",
        legend_style = {
            at = Coordinate(0.5, 1.05),
            anchor = "south",
            legend_columns = -1
        },
    },
    objects...
)

display(ax)
pgfsave("deflections.pdf", ax)
