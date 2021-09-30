# runtests.jl


for s in [
    "barrel_vault_examples",
    "clamped_hypar_examples", 
    "raasch_examples"
    ]
    include(s * ".jl")
end
