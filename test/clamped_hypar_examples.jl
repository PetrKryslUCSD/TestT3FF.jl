"""
From: CMES, vol.49, no.2, pp.81-110, 2009

The problem considered in this section is that of a hyperbolic paraboloid shell,
clamped along one side and free on three edges and loaded by self-weight
(Figure 16). This is a pure bending dominated problem and known to be a very
hard test for locking behaviour as suggested in References (Chapelle and Bathe,
1998; Bathe, Iosilevich, and Chapelle, 2000). 

Table 2 in Bathe, Iosilevich, and Chapelle, 2000 (computed with MITC16).
t/L   Strain energy  Displacement
1/100 1.6790e-3 9.3355e-5
1/1000 1.1013e-2 6.3941e-3
1/10,000 8.9867e-2 5.2988e-1

See also the table in  
An Improved Quadrilateral Flat Element with Drilling
Degrees of Freedom for Shell Structural Analysis
H. Nguyen-Van1 , N. Mai-Duy1 and T. Tran-Cong1
CMES, vol.49, no.2, pp.81-110, 2009.

The shell geometry is described by the
equation: z = x^2 −y^2 ; (x,y) ∈ − L/2 ; L/2
"""
module clamped_hypar_examples

using LinearAlgebra
using FinEtools
using FinEtoolsDeforLinear
using FinEtoolsFlexStructures.FESetShellT3Module: FESetShellT3
using FinEtoolsFlexStructures.FESetShellQ4Module: FESetShellQ4
using FinEtoolsFlexStructures.FEMMShellT3FFModule
using FinEtoolsFlexStructures.FEMMShellT3FFAModule
using FinEtoolsFlexStructures.RotUtilModule: initial_Rfield, linear_update_rotation_field!, update_rotation_field!
using FinEtoolsFlexStructures.VisUtilModule: plot_nodes, plot_midline, render, plot_space_box, plot_midsurface, space_aspectratio, save_to_json
using FinEtools.MeshExportModule.VTKWrite: vtkwrite

function _execute_full(tL_ratio = 1/100, g = 80*0.1^0, analyt_sol=-9.3355e-5, n = 32, visualize = false)
    # analytical solution for the vertical deflection and the midpoint of the
    # free edge 
    # Parameters:
    E = 2.0e11;
    nu = 0.3;
    L = 1.0
    thickness = tL_ratio * L
    # Bathe, Iosilevich, and Chapelle (2000) with a refined mesh of
    # high-order element MITC16
    formul = FEMMShellT3FFModule

    @info "Mesh: $n elements per side"

    tolerance = L/n/1000
    fens, fes = T3block(L,L,n,n);
    fens.xyz = xyz3(fens)
    for i in 1:count(fens)
        x=fens.xyz[i, 1]-L/2; y=fens.xyz[i, 2]-L/2;
        fens.xyz[i, :] .= (x, y, x^2 - y^2)
    end

    mater = MatDeforElastIso(DeforModelRed3D, E, nu)
    
    sfes = FESetShellT3()
    accepttodelegate(fes, sfes)
    femm = formul.make(IntegDomain(fes, TriRule(1), thickness), mater)
    # femm.drilling_stiffness_scale = 0.1
    # femm.mult_el_size = 0.2
    stiffness = formul.stiffness
    associategeometry! = formul.associategeometry!

    # Construct the requisite fields, geometry and displacement
    # Initialize configuration variables
    geom0 = NodalField(fens.xyz)
    u0 = NodalField(zeros(size(fens.xyz,1), 3))
    Rfield0 = initial_Rfield(fens)
    dchi = NodalField(zeros(size(fens.xyz,1), 6))

    # Apply EBC's
    # Clamped edge
    l1 = selectnode(fens; box = Float64[-L/2 -L/2 -Inf Inf -Inf Inf], inflate = tolerance)
    for i in [1,2,3,4,5,6]
        setebc!(dchi, l1, true, i)
    end
    applyebc!(dchi)
    numberdofs!(dchi);

    # Assemble the system matrix
    associategeometry!(femm, geom0)
    K = stiffness(femm, geom0, u0, Rfield0, dchi);

    # Midpoint of the free edge
    nl = selectnode(fens; box = Float64[L/2 L/2 0 0 -Inf Inf], inflate = tolerance)
    lfemm = FEMMBase(IntegDomain(fes, TriRule(1)))
    # computeforce!(forceout, XYZ, tangents, fe_label) = let
    #     n = cross(tangents[:, 1], tangents[:, 2])
    #     n = n / norm(n)
    #     forceout .= 0.0
    #     forceout[3] = -g * n[3]
    #     return forceout
    # end
    # fi = ForceIntensity(FFlt, 6, computeforce!)
    fi = ForceIntensity(FFlt[0, 0, -g, 0, 0, 0]);
    F = distribloads(lfemm, geom0, dchi, fi, 3);
    
    # Solve
    U = K\F
    scattersysvec!(dchi, U[:])
    targetu = dchi.values[nl, 3][1]
    @info "Solution: $(round(targetu/analyt_sol, digits = 4)*100)%"

        # Visualization
    if visualize
        scattersysvec!(dchi, (L/4)/abs(targetu).*U)
        update_rotation_field!(Rfield0, dchi)
        plots = cat(plot_space_box([[0 0 -L/2]; [L/2 L/2 L/2]]),
        #plot_nodes(fens),
            plot_midsurface(fens, fes; x = geom0.values, facecolor = "rgb(12, 12, 123)"),
            plot_midsurface(fens, fes; x = geom0.values, u = dchi.values[:, 1:3], R = Rfield0.values);
            dims = 1)
        pl = render(plots)
    end
    return targetu
end

function _execute_half(orientation = :a, tL_ratio = 1/100, g = 80*0.1^0, analyt_sol=-9.3355e-5, n = 32, visualize = false)
    # analytical solution for the vertical deflection and the midpoint of the
    # free edge 
    # Parameters:
    E = 2.0e11;
    nu = 0.3;
    L = 1.0
    thickness = tL_ratio * L
    # Bathe, Iosilevich, and Chapelle (2000) with a refined mesh of
    # high-order element MITC16
    formul = FEMMShellT3FFModule

    @info "Mesh: $n elements per side"

    tolerance = L/n/1000
    fens, fes = T3block(L,L/2,n,Int(round(n/2)),orientation);
    fens.xyz = xyz3(fens)
    for i in 1:count(fens)
        x=fens.xyz[i, 1]-L/2; y=fens.xyz[i, 2]-L/2;
        fens.xyz[i, :] .= (x, y, x^2 - y^2)
    end

    mater = MatDeforElastIso(DeforModelRed3D, E, nu)
    
    sfes = FESetShellT3()
    accepttodelegate(fes, sfes)
    femm = formul.make(IntegDomain(fes, TriRule(1), thickness), mater)
    # femm.drilling_stiffness_scale = 0.1
    # femm.mult_el_size = 0.2
    stiffness = formul.stiffness
    associategeometry! = formul.associategeometry!

    # Construct the requisite fields, geometry and displacement
    # Initialize configuration variables
    geom0 = NodalField(fens.xyz)
    u0 = NodalField(zeros(size(fens.xyz,1), 3))
    Rfield0 = initial_Rfield(fens)
    dchi = NodalField(zeros(size(fens.xyz,1), 6))

    # Apply EBC's
    # Clamped edge
    l1 = selectnode(fens; box = Float64[-L/2 -L/2 -Inf Inf -Inf Inf], inflate = tolerance)
    for i in [1,2,3,4,5,6]
        setebc!(dchi, l1, true, i)
    end
    # Symmetry Plane
    l1 = selectnode(fens; box = Float64[-Inf Inf 0 0 -Inf Inf], inflate = tolerance)
    for i in [2,4,6]
        setebc!(dchi, l1, true, i)
    end
    applyebc!(dchi)
    numberdofs!(dchi);

    # Assemble the system matrix
    associategeometry!(femm, geom0)
    K = stiffness(femm, geom0, u0, Rfield0, dchi);

    # Midpoint of the free edge
    nl = selectnode(fens; box = Float64[L/2 L/2 0 0 -Inf Inf], inflate = tolerance)
    lfemm = FEMMBase(IntegDomain(fes, TriRule(1)))
    fi = ForceIntensity(FFlt[0, 0, -g, 0, 0, 0]);
    F = distribloads(lfemm, geom0, dchi, fi, 3);
    
    # Solve
    U = K\F
    scattersysvec!(dchi, U[:])
    targetu = dchi.values[nl, 3][1]
    @info "Solution: $(round(targetu/analyt_sol, digits = 4)*100)%"

        # Visualization
    if visualize

        vtkwrite("clamped_hypar-$(orientation)-$(n).vtu", fens, fes; vectors = [("u", dchi.values[:, 1:3])])

        scattersysvec!(dchi, (L/4)/abs(targetu).*U)
        update_rotation_field!(Rfield0, dchi)
        plots = cat(plot_space_box([[0 0 -L/2]; [L/2 L/2 L/2]]),
        #plot_nodes(fens),
            plot_midsurface(fens, fes; x = geom0.values, facecolor = "rgb(12, 12, 123)"),
            plot_midsurface(fens, fes; x = geom0.values, u = dchi.values[:, 1:3], R = Rfield0.values);
            dims = 1)
        pl = render(plots)
    end
    return targetu
end

function test_convergence(orientation = :a)
    
    tL_ratios = [1/100, 1/1000, 1/10000]; 
    gs = [80*0.1^0, 80*0.1^1, 80*0.1^2]
    analyt_sols = [-9.3355e-5, -6.3941e-3, -5.2988e-1];
    
    ns = [4, 8, 16, 32, 64, 128, 256, 512]
    all_results = []
    for (tL_ratio, g, analyt_sol) in zip(tL_ratios, gs, analyt_sols)
        @info "Clamped hypar, t/L=$(tL_ratio)"
        results = Float64[]
        for n in ns
        # for n in [4, 8, 16, 32, 64, ]
            r = _execute_half(orientation, tL_ratio, g, analyt_sol, n, false)
            push!(results, r/analyt_sol)
        end   
        push!(all_results, results)
    end

    return ns, all_results
end

function test_0_001(orientation = :a, ns = [8, 16, 32, 48, ])
    tL_ratio = 1/1000 
    g = 80*0.1^1
    analyt_sol = -6.3941e-3
    
    
    @info "Clamped hypar, t/L=$(tL_ratio)"
    results = Float64[]
    for n in ns
        # for n in [4, 8, 16, 32, 64, ]
        r = _execute_half(orientation, tL_ratio, g, analyt_sol, n, true)
        push!(results, r/analyt_sol)
    end   
    

    return ns, results
end

# using Gnuplot

# ns = 1 ./ [4, 8, 16, 32, 64, 128, 256, ]
# results = [0.8839348674712617, 0.8750949157612452, 0.9199805802913757, 0.9619508790573108, 0.9856803572479892, 0.9955727622499687, 0.9993169031485688]   
# @gp ns results "with lp"      :-           
# results = [1.0429797613488236, 0.9314984628085947, 0.9365905225801154, 0.9565506281799385, 0.9764476699285441, 0.9902805329751646, 0.9968920296205528]   
# @gp :- ns results "with lp"    :-                             
# results = [1.251888013432877, 1.0155533090845452, 0.9678060658415124, 0.9718188010061173, 0.9812934066246979, 0.9889499817887738, 0.9953521405300628] 
# @gp :- ns results "with lp"



end # module

using PGFPlotsX

if true

using .clamped_hypar_examples

for orientation in (:a, :b)
    ns, all_results = clamped_hypar_examples.test_convergence(orientation)

    objects = []

    results = all_results[1]
    @pgf p = PGFPlotsX.Plot(
    {
    color = "black",
    line_width  = 0.7
    },
    Coordinates([v for v in  zip(1 ./ ns, results)])
    )
    push!(objects, p)
    push!(objects, LegendEntry("t/L=1/100"))

    results = all_results[2]
    @pgf p = PGFPlotsX.Plot(
    {
    color = "black",
    style = "dashed",
    line_width  = 0.7
    },
    Coordinates([v for v in  zip(1 ./ ns, results)])
    )
    push!(objects, p)
    push!(objects, LegendEntry("t/L=1/1000"))

    results = all_results[3]
    @pgf p = PGFPlotsX.Plot(
    {
    color = "black",
    style = "dotted",
    line_width  = 0.7
    },
    Coordinates([v for v in  zip(1 ./ ns, results)])
    )
    push!(objects, p)
    push!(objects, LegendEntry("t/L=1/10000"))

    @pgf ax = Axis(
    {
    xlabel = "Relative Element Size [ND]",
    ylabel = "Normalized Displacement [ND]",
        ymin = 0.7,
        ymax = 1.1,
    xmode = "log", 
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
    pgfsave("clamped_hypar_examples-dependence-on-t-L-$(orientation).pdf", ax)
end

end # if false

using PGFPlotsX

objects = []

# What is the source of these data points?
# Development of a strain-smoothed three-node triangular flat shell
# element with drilling degrees of freedom
# C.M. Shin, B.C. Lee / Finite Elements in Analysis and Design 86 (2014) 71–80
# For some reason the results there are scaled with -5.958e-3.


# Allman = [0.07, 0.3074, 0.829, 0.9545] .* (-5.958e-3) ./ (-6.3941e-3)
# Providas_Kattis = [1.07, 1.00, 0.9859, 0.9857] .* (-5.958e-3) ./ (-6.3941e-3)
# Cook_flat_stiffened = [0.5632, 0.9399, 0.9849, 0.9913] .* (-5.958e-3) ./ (-6.3941e-3)
# Shin_Lee = [1.0826 1.0144 0.9986 0.9991]  .* (-5.958e-3) ./ (-6.3941e-3)

# all_results = [("Present", ns, results, "*"), ("Allman", Allman, "x"), ("Providas, Kattis", Providas_Kattis, "triangle"), ("Cook, flat, stiffened", Cook_flat_stiffened, "square"), ("Shin, Lee", Shin_Lee, "diamond")]

# Here we take the reference results from Bathe, Iosilevich, Chapelle.

# An efficient three‑node triangular Mindlin–Reissner flat shell element
# Hosein Sangtarash1 · Hamed Ghohani Arab1 · Mohammad R. Sohrabi1 · Mohammad R. Ghasemi1
# TMRFS = ("TMRFS", [4, 8, 16, 24], [0.533 0.857 0.990 0.998], "diamond")
# MITC3p = ("MITC3+", [4, 8, 16, ], [0.407 0.768 0.93], "triangle")
# all_results = [("Present", ns, results, "*"), TMRFS, MITC3p]

# Performance of the MITC3+ and MITC4+ shell elements in widely-used
# benchmark problems
# Yeongbin Ko a, Youngyu Lee b, Phill-Seung Lee a,⇑, Klaus-Jürgen Bathe c
# a Department of Mechanical Engineering, Korea Advanced Institute of Science and Technology, 291 Daehak-ro, Yuseong-gu, 
MITC3p = ("MITC3+", [4, 8, 16, 32, 64], [0.9533 0.9589 0.9728 0.9868 0.9951], "diamond")

@show ns, resultsa = clamped_hypar_examples.test_0_001(:a, [4, 8, 16, 32, 64, 128, 256])
@show ns, resultsb = clamped_hypar_examples.test_0_001(:b, [4, 8, 16, 32, 64, 128, 256])
all_results = [("Present (a)", ns, resultsa, "*"), ("Present (b)", ns, resultsb, "o"), MITC3p]


for r in  all_results
    @pgf p = PGFPlotsX.Plot(
    {
    color = "black",
    line_width  = 0.7, 
    style = "solid",
    mark = "$(r[4])"
    },
    Coordinates([v for v in  zip(1 ./ (r[2]), abs.(1 .- r[3]))])
    )
    push!(objects, p)
    push!(objects, LegendEntry("$(r[1])"))
end


@pgf ax = Axis(
    {
        xlabel = "Relative element size [ND]",
        ylabel = "Normalized Displacement Error [ND]",
        # xmin = range[1],
        # xmax = range[2],
        xmode = "log", 
        ymode = "log",
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
pgfsave("clamped_hypar_examples-0_001-errors.pdf", ax)


using .clamped_hypar_examples
ns, all_results = clamped_hypar_examples.test_0_001(:a, 4 .* [24, 48, 96])
q1, q2, q3 = all_results[:]
qtrue = (q2^2 - q1 * q3) / (2*q2 - q1 - q3)
@show (q1, q2, q3) .* -6.3941, qtrue .* -6.3941
ns, all_results = clamped_hypar_examples.test_0_001(:b, 4 .* [24, 48, 96])
q1, q2, q3 = all_results[:]
qtrue = (q2^2 - q1 * q3) / (2*q2 - q1 - q3)
@show (q1, q2, q3) .* -6.3941, qtrue .* -6.3941

q1, q2, q3 = MITC3p[3][3:end]
qtrue = (q2^2 - q1 * q3) / (2*q2 - q1 - q3)
@show (q1, q2, q3) .* -6.3941, qtrue .* -6.3941
