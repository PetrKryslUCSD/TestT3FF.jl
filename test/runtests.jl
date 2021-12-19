# runtests.jl


for s in readdir("."; join=true, sort=true)
    f, e = splitext(s)
    if e == ".pdf" || e == ".vtu"
  try
   rm(s)
        catch
            @warn "Failed to remove $s"
        end
    end
end

for s in [
    "LE5_Z_cantilever_examples",
    "convergence_metaanalysis",
    "cos_2t_press_hyperboloid_free_examples",
    "vis_scordelis_lo_examples",
    "barrel_vault_examples",
    "cos_2t_press_cylinder_fixed_examples",
    "hemisphere_open_examples",
    "scordelis_lo_examples",
    "barrel_w_stiffeners_examples",
    "cos_2t_press_cylinder_free_examples",
    "pinched_cylinder_examples",
    "clamped_hypar_examples",
    "cos_2t_press_hyperboloid_fixed_examples",
    "raasch_examples",
    "ss_circular_plate_udl_examples",
    ]
    include(s * ".jl")
end
