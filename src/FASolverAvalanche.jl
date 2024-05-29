module FASolverAvalanche

include("./HelperFunctions/meshParser.jl") # For Parsing input files 
include("./HelperFunctions/mesher.jl") # Mesh Computations
include("./HelperFunctions/preProcess.jl") # Pre-processing and structuring mesh 
include("./initialConditions.jl") # Initialize discretized system 
include("./solver.jl") # Solve the discretized system 


end
