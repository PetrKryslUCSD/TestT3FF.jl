# TestT3FF.jl

Illustration of the T3FF general flat-facet shell element capabilities.
The paper describing this element has been submitted 12/2021.

Run
```
include("top.jl")
```
and then execute
```
using Pkg; Pkg.test(); 
```
The `test` directory will have generated files, `.vtu` files for visualization with Paraview, and `.pdf` files for the generated graphs.
Beware: running the above will first clean all the output files from the previous run.