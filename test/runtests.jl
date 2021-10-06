# runtests.jl


for s in [
    "barrel_vault_examples",
    "clamped_hypar_examples", 
    "pinched_cylinder_examples",
    "ss_circular_plate_udl_examples",
    "hemisphere_open_examples", 
    "raasch_examples"
    ]
    include(s * ".jl")
end
