# TestT3FF.jl

Illustration of the T3FF general flat-facet shell element capabilities.
The paper describing this element has been submitted 12/2021.

The functionality of the shell finite element is provided by the package
[`FinEtoolsFlexStructures`]
(https://github.com/PetrKryslUCSD/FinEtoolsFlexStructures.jl).

Assuming the `TestT3FF.jl` package was cloned from Github, change your working directory
into the `TestT3FF.jl` folder, start `julia`, and then run
```
include("top.jl")
```
To produce the results presented in the paper, execute
```
using Pkg; Pkg.test(); 
```
The `test` directory will hold the generated files, `.vtu` files for 
visualization with Paraview, and `.pdf` files for the generated graphs.
Beware: running the above will first clean all the output files from the previous run.