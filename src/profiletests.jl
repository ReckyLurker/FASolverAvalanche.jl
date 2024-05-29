using Profile

include("./Mesher.jl")
include("./MeshParser.jl")
include("./preprocess.jl")

# Mesh Generation and text Parsing

# Compilation
mesh, points, faces = generateMesh("/home/recklurker/RWTHIntern/Julia/points",
                    "/home/recklurker/RWTHIntern/Julia/faces",
                    "/home/recklurker/RWTHIntern/Julia/faceLabels", visualizeMesh=true)

# Pure runtime
@time for i = 1:10
    mesh, points, faces = generateMesh("/home/recklurker/RWTHIntern/Julia/points",
                    "/home/recklurker/RWTHIntern/Julia/faces",
                    "/home/recklurker/RWTHIntern/Julia/faceLabels", visualizeMesh=true)
end