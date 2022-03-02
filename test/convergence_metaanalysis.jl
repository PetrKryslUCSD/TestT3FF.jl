results =  []
which = 0
mypush!(a) = let
    global which
    which += 1
    push!(results, (which, a))
end

# Shin, Lee
# mypush!([100.46, 98.71, 98.74])
# mypush!([90.75, 92.90, 95.06])
mypush!([73.45, 81.47, 87.36])
mypush!([108.3, 103.04, 101.39])
mypush!([101.85, 100.79, 100.43])
mypush!([102.31, 100.84, 100.43])

# Sangtarash et al
mypush!([104.5, 101.0, 100.3])
# mypush!([0.94, 0.98, 1.01] .* 100)
mypush!([104.8, 100.5, 99.7])
mypush!([94.3, 97.2, 98.8])
mypush!([103.6, 100.2, 99.6])
mypush!([104.7, 100.5, 99.7])
mypush!([116.8, 102.8, 100.8])
# mypush!([77.1, 83.0, 94.0])
mypush!([94.3, 97.2, 98.9])
mypush!([106.6, 104.4, 103.4])
# mypush!([99.5, 98.6, 99.3])
mypush!([93.7, 97.4, 99.0])
mypush!([94.3, 97.3, 98.9])
# mypush!([])
mypush!([99.4, 99.7, 99.9])

# Rezaiee
mypush!([58.08, 90.94, 98.84])
mypush!([52.86, 88.31, 98.09])
mypush!([44.63, 84.28, 96.68])
mypush!([62.04, 92.25, 99.04])
# mypush!([# ])
mypush!([94.37, 97.25, 98.90])
mypush!([104.7, 101.1, 99.6])

# Briassoulis
mypush!([0.964, 0.984, 0.998] .*100)
mypush!([120.1, 104.6, 101.0])
# mypush!([94.0, 97.0, 100.0]
# mypush!([])
mypush!([108.3, 101.5, 100.0])
# mypush!([94.0, 98.0, 101.0])
mypush!([95.0, 97.4, 98.7])

# Argyris
# mypush!([69.7, 90.2, 100.1])

# Cen, LI
mypush!([103.74, 100.29, 99.61])
mypush!([100.14, 99.56, 99.51])
mypush!([103.56, 100.19, 99.66])
mypush!([116.8, 102.8, 100.8])
# mypush!([92.84, 96.09, 99.08])
# mypush!([77.12, 83.07, 94.31])
mypush!([119.12, 104.2, 100.63])
# mypush!([])
# mypush!([94.2, 100.8, 100.5])


# Ko et al
mypush!([95.50, 98.51, 99.32])
mypush!([94.32, 97.26, 98.86])
mypush!([104.8, 100.5, 99.73])
mypush!([108.3, 101.5, 100])
mypush!([104.8, 100.5, 99.6])
mypush!([104.7, 100.5, 99.74])

# Jun, et al
# mypush!([73.12, 87.43, 95.93])
# mypush!([66.77, 85.58, 95.40])
# mypush!([72.78, 90.15, 100.43])
mypush!([96.10, 99.31, 99.83])
mypush!([89.22, 97.62, 99.50])
mypush!([94.31, 97.29, 98.88])
# mypush!([94.01, 98.02, 100.99])
# mypush!([104.99, 100.78])
mypush!([121.89, 105.39, 101.69])

# Nguyen-Van et al
mypush!([119.12, 104.2, 100.63])
mypush!([119.25, 104.22, 100.66])
mypush!([120.1, 104.6, 101.0])
mypush!([97.6, 98.6, 99.3])
mypush!([104.7, 100.5, 99.7])
mypush!([121.9, 105.4, 101.7])

# Mostafa
mypush!([104.9, 101.3, 99.6])
mypush!([104.7, 101.1, 99.6])
# mypush!([99.7, ])
mypush!([98.6, 99.3, 99.6])
mypush!([96.0, 98.4, 99.9])
mypush!([102.9, 100.1, 99.2])

# Coupling Effect
mypush!([.326802, .316916, .313979] ./ (0.3024) .* 100)
mypush!([0.329038, 0.317741, 0.314225] ./ (0.3024) .* 100)
# # mypush!([0.095964, ] ./ (0.3024) .* 100)
# mypush!([.220731, .274365, .307401] ./ (0.3024) .* 100)
# mypush!([.284311, .293325, .305336] ./ (0.3024) .* 100)
mypush!([.316522, .303931, .301618] ./ (0.3024) .* 100)

mypush!(-0.23851E+00, -2.994369E-01, -3.021210E-01] ./ (-0.3024) .* 100)

# The extrapolated value with Richardson formula is 0.3022447

differences = []
for r in results
    @show r[1]
    @show q1, q2, q3 = r[2]
    @show qtrue = (q2^2 - q1 * q3) / (2*q2 - q1 - q3)
    push!(differences, (r[1], qtrue - 100))
end

differences = sort(differences, lt = (x, y) -> x[2] < y[2])
@show differences
# @show [(i, v) in  for (i, v) in enumerate(differences)]

using PGFPlotsX

objects = []


@pgf p = PGFPlotsX.Plot(
{
color = "black",
line_width  = 0.7, 
style = "solid",
mark = "*"
},
Coordinates([(i, v[2]) for (i, v) in enumerate(differences)])
)
push!(objects, p)
# mypush!(objects, LegendEntry("$(r[1])"))


@pgf ax = Axis(
{
xlabel = "Study [ND]",
ylabel = "Departure from 0.3024 [\\%]",
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
pgfsave("Convergence-metaanalysis.pdf", ax)
