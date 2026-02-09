"""
Pinched cylinder with diagphram supports and concentrated force
"""
module pinched_cylinder_examples

using FinEtools
using FinEtoolsDeforLinear
using FinEtoolsFlexStructures.FESetShellT3Module: FESetShellT3
using FinEtoolsFlexStructures.FEMMShellT3FFModule
using FinEtoolsFlexStructures.RotUtilModule: initial_Rfield, update_rotation_field!
using VisualStructures: plot_nodes, plot_midline, render, plot_space_box, plot_midsurface, space_aspectratio, save_to_json

function _execute(n = 2, visualize = true)
    E = 3e6;
    nu = 0.3;
    thickness = 3.0;
    # analytical solution for the vertical deflection under the load
    analyt_sol=-1.82488e-5;
    formul = FEMMShellT3FFModule

    @info "Mesh: $n elements per side"

    # Mesh
    R = 300.0;
    L = 600.0;
    
    tolerance = R/n/1000
    fens, fes = T3block(90/360*2*pi,L/2,n,n);
    fens.xyz = xyz3(fens)
    for i in 1:count(fens)
        a=fens.xyz[i, 1]; y=fens.xyz[i, 2];
        fens.xyz[i, :] .= (R*sin(a), y, R*cos(a))
    end
    
    mater = MatDeforElastIso(DeforModelRed3D, E, nu)

    
    sfes = FESetShellT3()
    accepttodelegate(fes, sfes)
    femm = formul.make(IntegDomain(fes, TriRule(1), thickness), mater)
    femm.drilling_stiffness_scale = 0.1
    associategeometry! = formul.associategeometry!
    stiffness = formul.stiffness

    # Construct the requisite fields, geometry and displacement
    # Initialize configuration variables
    geom0 = NodalField(fens.xyz)
    u0 = NodalField(zeros(size(fens.xyz,1), 3))
    Rfield0 = initial_Rfield(fens)
    dchi = NodalField(zeros(size(fens.xyz,1), 6))

    # Apply EBC's
    # rigid diaphragm
    l1 = selectnode(fens; box = Float64[-Inf Inf 0 0 -Inf Inf], inflate = tolerance)
    for i in [1,3]
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
    # plane of symmetry perpendicular to Z
    l1 = selectnode(fens; box = Float64[-Inf Inf -Inf Inf 0.0 0.0], inflate = tolerance)
    for i in [3,4,5]
        setebc!(dchi, l1, true, i)
    end
    applyebc!(dchi)
    numberdofs!(dchi);

    # Assemble the system matrix
    associategeometry!(femm, geom0)
    K = stiffness(femm, geom0, u0, Rfield0, dchi);

    # Load
    nlc = selectnode(fens; box = Float64[0 0 L/2 L/2 -Inf Inf], inflate = tolerance)
    # If desired, increase the tolerance, for instance to half of the thickness.
    nl = selectnode(fens; box = Float64[0 0 L/2 L/2 -Inf Inf], inflate = tolerance)
    nln = length(nl)
    loadbdry = FESetP1(reshape(nl, nln, 1))
    lfemm = FEMMBase(IntegDomain(loadbdry, PointRule()))
    fi = ForceIntensity(FFlt[0, 0, -1/4/nln, 0, 0, 0]);
    F = distribloads(lfemm, geom0, dchi, fi, 3);

    # Solve
    solve_blocked!(dchi, K, F)
    U = gathersysvec(dchi, DOF_KIND_ALL)
    @show n,  dchi.values[nlc, 3][1]/analyt_sol*100

    # Visualization
    if visualize
        scattersysvec!(dchi, (L/8)/maximum(abs.(U)).*U)
        update_rotation_field!(Rfield0, dchi)
        plots = cat(plot_space_box([[0 0 -L/2]; [L/2 L/2 L/2]]),
        #plot_nodes(fens),
            plot_midsurface(fens, fes; x = geom0.values, u = dchi.values[:, 1:3], R = Rfield0.values);
            dims = 1)
        pl = render(plots)
    end
    return dchi.values[nlc, 3][1]/analyt_sol*100
end

function test_convergence()
    @info "Pinched cylinder"
    ns = [4, 8, 16, 24, 32] 
    results = Float64[]
    for n in ns 
        v = _execute(n, false)
        push!(results, v)
    end
    return ns, results
end

end # module

using .pinched_cylinder_examples
ns, results = pinched_cylinder_examples.test_convergence()

# [ Info: Mesh: 4 elements per side                                          
# (n, (dchi.values[nl, 3] / analyt_sol) * 100) = (4, [41.55321964302738])    
# [ Info: Mesh: 8 elements per side                                          
# (n, (dchi.values[nl, 3] / analyt_sol) * 100) = (8, [76.87858207030617])    
# [ Info: Mesh: 12 elements per side                                         
# (n, (dchi.values[nl, 3] / analyt_sol) * 100) = (12, [87.73634899399313])   
# [ Info: Mesh: 16 elements per side                                         
# (n, (dchi.values[nl, 3] / analyt_sol) * 100) = (16, [92.41804796585181])   
# [ Info: Mesh: 24 elements per side                                         
# (n, (dchi.values[nl, 3] / analyt_sol) * 100) = (24, [96.49232304717079])   
# [ Info: Mesh: 48 elements per side                                         
# (n, (dchi.values[nl, 3] / analyt_sol) * 100) = (48, [99.6766661386644])    
# [ Info: Mesh: 96 elements per side                                         
# (n, (dchi.values[nl, 3] / analyt_sol) * 100) = (96, [100.81626349991888])  
# [ Info: Mesh: 128 elements per side                                        
# (n, (dchi.values[nl, 3] / analyt_sol) * 100) = (128, [101.05826709146844]) 
# [ Info: Mesh: 256 elements per side                                        
# (n, (dchi.values[nl, 3] / analyt_sol) * 100) = (256, [101.43102704824867]) 
# [ Info: Mesh: 400 elements per side                                        
# (n, (dchi.values[nl, 3] / analyt_sol) * 100) = (400, [101.6018373070541])  
# [ Info: Mesh: 600 elements per side                                        
# (n, (dchi.values[nl, 3] / analyt_sol) * 100) = (600, [101.73923429243655]) 

# These results come from Shin, Lee (2014), Table 3.

# Mesh Subdivision
# 4, 8, 16, 24, 32
Allman = [
59.02
92.41
100.43
100.53
100.29
] .* (1.8541/1.82488)
Cook_flat_stiffened = [
53.71
89.73
99.67
100.19
100.10
] .* (1.8541/1.82488)
Providas_Kattis = [
45.35
85.64
98.24
99.52
99.74
] .* (1.8541/1.82488)

using PGFPlotsX

objects = []



all_results = [("Present", results, "*"), ("Allman", Allman, "x"), ("Providas, Kattis", Providas_Kattis, "triangle"), ("Cook, flat, stiffened", Cook_flat_stiffened, "square")]

for r in  all_results
    @pgf p = PGFPlotsX.Plot(
    {
    color = "black",
    line_width  = 0.7, 
    style = "solid",
    mark = "$(r[3])"
    },
    Coordinates([v for v in  zip(ns, r[2])])
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
        legend_style = {
            at = Coordinate(0.5, 1.05),
            anchor = "south",
            legend_columns = -1
        },
    },
    objects...
)

display(ax)
pgfsave("pinched_cylinder-convergence.pdf", ax)
