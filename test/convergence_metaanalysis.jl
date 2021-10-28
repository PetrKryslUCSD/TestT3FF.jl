results =  []

# Shin, Lee
# push!(results, [100.46, 98.71, 98.74])
# push!(results, [90.75, 92.90, 95.06])
push!(results, [73.45, 81.47, 87.36])
push!(results, [108.3, 103.04, 101.39])
push!(results, [101.85, 100.79, 100.43])
push!(results, [102.31, 100.84, 100.43])

# Sangtarash et al
push!(results, [104.5, 101.0, 100.3])
# push!(results, [0.94, 0.98, 1.01] .* 100)
push!(results, [104.8, 100.5, 99.7])
push!(results, [94.3, 97.2, 98.8])
push!(results, [103.6, 100.2, 99.6])
push!(results, [104.7, 100.5, 99.7])
push!(results, [116.8, 102.8, 100.8])
# push!(results, [77.1, 83.0, 94.0])
push!(results, [94.3, 97.2, 98.9])
push!(results, [106.6, 104.4, 103.4])
# push!(results, [99.5, 98.6, 99.3])
push!(results, [93.7, 97.4, 99.0])
push!(results, [94.3, 97.3, 98.9])
# push!(results, [])
push!(results, [99.4, 99.7, 99.9])

# Rezaiee
push!(results, [58.08, 90.94, 98.84])
push!(results, [52.86, 88.31, 98.09])
push!(results, [44.63, 84.28, 96.68])
push!(results, [62.04, 92.25, 99.04])
# push!(results, [# ])
push!(results, [94.37, 97.25, 98.90])
push!(results, [104.7, 101.1, 99.6])

# Briassoulis
push!(results, [0.964, 0.984, 0.998] .*100)
push!(results, [120.1, 104.6, 101.0])
# push!(results, [94.0, 97.0, 100.0])
push!(results, [104.8, 100.5, 99.6])
# push!(results, [])
push!(results, [108.3, 101.5, 100.0])
# push!(results, [94.0, 98.0, 101.0])
push!(results, [95.0, 97.4, 98.7])

# Argyris
# push!(results, [69.7, 90.2, 100.1])

# Cen, LI
push!(results, [103.74, 100.29, 99.61])
push!(results, [100.14, 99.56, 99.51])
push!(results, [103.56, 100.19, 99.66])
push!(results, [116.8, 102.8, 100.8])
# push!(results, [92.84, 96.09, 99.08])
# push!(results, [77.12, 83.07, 94.31])
push!(results, [119.12, 104.2, 100.63])
# push!(results, [])
# push!(results, [94.2, 100.8, 100.5])


# Ko et al
push!(results, [95.50, 98.51, 99.32])
push!(results, [94.32, 97.26, 98.86])
push!(results, [104.8, 100.5, 99.73])
push!(results, [108.3, 101.5, 100])
push!(results, [104.8, 100.5, 99.6])
push!(results, [104.7, 100.5, 99.74])

# Jun, et al
push!(results, [73.12, 87.43, 95.93])
push!(results, [66.77, 85.58, 95.40])
# push!(results, [72.78, 90.15, 100.43])
push!(results, [96.10, 99.31, 99.83])
push!(results, [89.22, 97.62, 99.50])
push!(results, [94.31, 97.29, 98.88])
# push!(results, [94.01, 98.02, 100.99])
# push!(results, [104.99, 100.78])
push!(results, [121.89, 105.39, 101.69])


differences = []
for r in results
    q1, q2, q3 = r
    qtrue = (q2^2 - q1 * q3) / (2*q2 - q1 - q3)
    push!(differences, qtrue - 100)
end

differences = sort(differences)

using PGFPlotsX

objects = []


@pgf p = PGFPlotsX.Plot(
{
color = "black",
line_width  = 0.7, 
style = "solid",
mark = "x"
},
Coordinates([v for v in enumerate(differences)])
)
push!(objects, p)
# push!(objects, LegendEntry("$(r[1])"))


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
