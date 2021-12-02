"""
The barrel vault (Scordelis-Lo) roof is one of the benchmarks for linear elastic
analysis of shells. 

The candidate element's usefulness in irregular geometries (and most practical
cases involve a high degree of geometric irregularity) is tested. As would be
expected,the irregular mesh results are not as good as those provided by a
regular meshwith the same number of variables. 

Problem description

The physical basis of the problem is a deeply arched roof supported only
bydiaphragms at its curved edges (an aircraft hanger), deforming under its own
weight. It is interesting to observe that the geometry is such that the
centerpoint of the roof moves upward under the self-weight(downwardly directed)
load. Perhaps this is one reason why the problem is not straightforward
numerically. 
"""
module scordelis_lo_examples

using LinearAlgebra
using FinEtools
using FinEtoolsDeforLinear
using FinEtoolsFlexStructures.FESetShellT3Module: FESetShellT3
using FinEtoolsFlexStructures.FESetShellQ4Module: FESetShellQ4
using FinEtoolsFlexStructures.FEMMShellT3FFModule
using FinEtoolsFlexStructures.RotUtilModule: initial_Rfield, linear_update_rotation_field!, update_rotation_field!
using FinEtoolsFlexStructures.VisUtilModule: plot_nodes, plot_midline, render, plot_space_box, plot_midsurface, space_aspectratio, save_to_json

using Infiltrator

function _execute_dsg_model(formul, n = 8, visualize = true)
    # analytical solution for the vertical deflection and the midpoint of the
    # free edge 
    analyt_sol=-0.3024;
    # Parameters:
    E=4.32e8;
    nu=0.0;
    thickness = 0.25; # geometrical dimensions are in feet
    R = 25.0;
    L = 50.0;
    
    tolerance = R/n/1000
    fens, fes = T3block(40/360*2*pi,L/2,n,n);
    fens.xyz = xyz3(fens)
    for i in 1:count(fens)
        a=fens.xyz[i, 1]; y=fens.xyz[i, 2];
        fens.xyz[i, :] .= (R*sin(a), y, R*(cos(a)-1))
    end

    mater = MatDeforElastIso(DeforModelRed3D, E, nu)
    
    sfes = FESetShellT3()
    accepttodelegate(fes, sfes)
    femm = formul.make(IntegDomain(fes, TriRule(1), thickness), mater)
    stiffness = formul.stiffness
    associategeometry! = formul.associategeometry!

    # Construct the requisite fields, geometry and displacement
    # Initialize configuration variables
    geom0 = NodalField(fens.xyz)
    u0 = NodalField(zeros(size(fens.xyz,1), 3))
    Rfield0 = initial_Rfield(fens)
    dchi = NodalField(zeros(size(fens.xyz,1), 6))

    # Apply EBC's
    # rigid diaphragm
    l1 = selectnode(fens; box = Float64[-Inf Inf 0 0 -Inf Inf], inflate = tolerance)
    for i in [1,3,5]
        setebc!(dchi, l1, true, i)
    end
    # plane of symmetry perpendicular to Y
    l1 = selectnode(fens; box = Float64[-Inf Inf L/2 L/2 -Inf Inf], inflate = tolerance)
    for i in [2,4,6]
        setebc!(dchi, l1, true, i)
    end
    # plane of symmetry perpendicular to X
    l1 = selectnode(fens; box = Float64[0 0 -Inf Inf -Inf Inf], inflate = tolerance)
    for i in [1,5,6]
        setebc!(dchi, l1, true, i)
    end
    applyebc!(dchi)
    numberdofs!(dchi);

    # Assemble the system matrix
    associategeometry!(femm, geom0)
    K = stiffness(femm, geom0, u0, Rfield0, dchi);

    # Midpoint of the free edge
    nl = selectnode(fens; box = Float64[sin(40/360*2*pi)*25 sin(40/360*2*pi)*25 L/2 L/2 -Inf Inf], inflate = tolerance)
    lfemm = FEMMBase(IntegDomain(fes, TriRule(3)))
    fi = ForceIntensity(FFlt[0, 0, -90, 0, 0, 0]);
    F = distribloads(lfemm, geom0, dchi, fi, 3);
    
    # Solve
    U = K\F
    scattersysvec!(dchi, U[:])
    result =   dchi.values[nl, 3][1]
    @info "Solution: $(result), $(round(result/analyt_sol*100, digits = 4))%"

    # Visualization
    if visualize
        scattersysvec!(dchi, (L/8)/maximum(abs.(U)).*U)
        update_rotation_field!(Rfield0, dchi)
        plots = cat(plot_space_box([[0 0 -L/2]; [L/2 L/2 L/2]]),
            #plot_nodes(fens),
            plot_midsurface(fens, fes; x = geom0.values, facecolor = "rgb(12, 12, 123)"),
            plot_midsurface(fens, fes; x = geom0.values, u = dchi.values[:, 1:3], R = Rfield0.values);
            dims = 1)
        pl = render(plots)
    end

    result
end

function test_convergence(ns = [4, 6, 8, 16, 32])
    formul = FEMMShellT3FFModule
    @info "Scordelis-Lo shell"
    results = []
    for n in ns
        v = _execute_dsg_model(formul, n, false)
        push!(results, v/(-0.3024)*100)
    end
    return ns, results
end

end # module

using .scordelis_lo_examples

ns, results = scordelis_lo_examples.test_convergence()


# These results come from Table 9 of An efficient three‑node triangular
# Mindlin–Reissner flat shell element, Hosein Sangtarash1 · Hamed Ghohani
# Arab1 · Mohammad R. Sohrabi1 · Mohammad R. Ghasemi1
# they are all normalized relative to 0.3024# .   # 
# Mesh Subdivision
# 4, 8, 16, 32

Allman = [
1.004
0.987
0.987
0.988] .* 100
Cook = [
0.907
0.929
0.950
0.981] .* 100
Providas_Kattis = [
0.734
0.815
0.873
0.967] .* 100
Shin_Lee = [
1.379
1.023
1.004
missing] .* 100
MITC3plus = [
0.669
missing
0.857
0.955] .* 100
TMRFS = [
0.924
0.963
0.974
0.998
] .* 100

using PGFPlotsX

objects = []

ns, results = scordelis_lo_examples.test_convergence()

compensate(r) = (0.3024/0.3006) .* r
all_results = [("Present", compensate(results), "*"), ("Allman", compensate(Allman), "x"), ("ProvKat", compensate(Providas_Kattis), "triangle"), ("Cook", compensate(Cook), "square"), ("SL", compensate(Shin_Lee), "o"), ("MITC3+", compensate(MITC3plus), "diamond")]

for r in  all_results
    @pgf p = PGFPlotsX.Plot(
    {
    color = "black",
    line_width  = 0.7, 
    style = "solid",
    mark = "$(r[3])"
    },
    Coordinates([v for v in  zip(ns, r[2]) if v[2] !== missing])
    )
    push!(objects, p)
    push!(objects, LegendEntry("$(r[1])"))
end


@pgf ax = Axis(
    {
        xlabel = "Number of Elements / side [ND]",
        ylabel = "Normalized Displacement [ND]",
        # xmin = range[1],
        # xmax = range[2],
        xmode = "linear", 
        ymode = "linear",
        yminorgrids = "true",
        grid = "both",
        legend_pos = "north east"
    },
    objects...
)

display(ax)
pgfsave("scordelis_lo_examples-convergence.pdf", ax)

# ns, results = scordelis_lo_examples.test_convergence([4, 8, 16])
# q1, q2, q3 = results
# @show qtrue = (q2^2 - q1 * q3) / (2*q2 - q1 - q3)

# ns, results = scordelis_lo_examples.test_convergence(2 .* [4, 8, 16])
# q1, q2, q3 = results
# @show qtrue = (q2^2 - q1 * q3) / (2*q2 - q1 - q3)

ns, results = scordelis_lo_examples.test_convergence([16, 32, 64])
q1, q2, q3 = results
@show qtrue = (q2^2 - q1 * q3) / (2*q2 - q1 - q3)  *  0.3024/100

# ns, results = scordelis_lo_examples.test_convergence(8 .* [4, 8, 16])
# q1, q2, q3 = results
# @show qtrue = (q2^2 - q1 * q3) / (2*q2 - q1 - q3)

# ns, results = scordelis_lo_examples.test_convergence(16 .* [4, 8, 16])
# q1, q2, q3 = results
# @show qtrue = (q2^2 - q1 * q3) / (2*q2 - q1 - q3)

# ns, results = scordelis_lo_examples.test_convergence(32 .* [4, 8, 16])
# q1, q2, q3 = results
# @show qtrue = (q2^2 - q1 * q3) / (2*q2 - q1 - q3)

# ns, results = scordelis_lo_examples.test_convergence(64 .* [4, 8, 16])
# q1, q2, q3 = results
# @show qtrue = (q2^2 - q1 * q3) / (2*q2 - q1 - q3)

# julia> q1, q2, q3 = results                                                             
# 3-element Vector{Any}:                                                                  
#  99.65925843568742                                                                      
#  99.6820694134264                                                                       
#  99.68801945695984                                                                      
                                                                                        
# julia> @show qtrue = (q2^2 - q1 * q3) / (2*q2 - q1 - q3)                                
# qtrue = (q2 ^ 2 - q1 * q3) / ((2q2 - q1) - q3) = 99.69011916380857                      
# 99.69011916380857    

# using FinEtools.AlgoBaseModule 
# AlgoBaseModule.richextrapol(Float64.(results) , 1 ./ ns)                         
# (99.69011916375054, 1.9387561566847429, 1440.1118605939823, [3.750472155061857e-14, 3.239596091386687e-14, -1.0034074265918846e-14])  

# # Data from Kiendl
# # Number of knots
# #    4                   8                   16                  32                  64
# 2  0.069489140101432   0.238098983889902   0.295702595727846   0.300237752650471   0.300558363037311
# 3  0.279944686988951   0.300065498003107   0.300584157108943   0.300592331224518   0.300592454396341
# 4  0.300000190520180   0.300589846259349   0.300592432120400   0.300592456541168   0.300592456632799
# 5  0.300966825455877   0.300589733200638   0.300592456543583   0.300592456635101   0.300592456617171

# # Suitable for extrapolation: p=5, first three results
# # Surprising: convergence from above?
# q1, q2, q3 = 0.300966825455877,   0.300589733200638,   0.300592456543583 
# @show qtrue = (q2^2 - q1 * q3) / (2*q2 - q1 - q3)  # = 0.30059243701673916 
# # Not suitable for extrapolation: p=5, last three results
# q1, q2, q3 = 0.300592456543583,   0.300592456635101,   0.300592456617171
# @show qtrue = (q2^2 - q1 * q3) / (2*q2 - q1 - q3)