include("./MeshParser.jl")
include("./Mesher.jl")
include("./preprocess.jl")
import GLMakie as Mke # Visualizer
import Images.RGB as RGB # Color Schemes
import Images.Gray as Gray 

# Linear Map from range(a,b) to range(c,d) [Assuming x lies in (a,b)]
function map(a,b,c,d,x)
    return c + (d - c)/(b-a)*(x-a)
end

mesh, points, faces = generateMesh("/home/recklurker/RWTHIntern/Julia/points",
                    "/home/recklurker/RWTHIntern/Julia/faces",
                    "/home/recklurker/RWTHIntern/Julia/faceLabels", visualizeMesh=true)

# To Visualize mesh
# viz(mesh, color=:white, showsegments=true)

# Calculate Normals
# normals = calculateNormals(points, faces)
# center_normals = [p[1] for p in normals] # Center of Face

# areas = computeFaceAreas(points, faces, center_normals)
# println(areas)
# direction_normals = [p[2] for p in normals] # Direction of Normals

# Visualize Normals as line segments
# normal_endpoints = [center_normals[i] + 10*direction_normals[i] for i in 1:length(normals)]
# line_normals = [Segment(center_normals[i], normal_endpoints[i]) for i in 1:length(normals)]
# viz(line_normals)

# Visualize Normals as Colormaps
# color_normals = [RGB(map(-1,1,0,1,direction_normals[i].coords[1]), 
#                  map(-1,1,0,1,direction_normals[i].coords[2]), 
#                  map(-1,1,0,1,direction_normals[i].coords[3])) for i in 1:length(normals)]

# viz(mesh, color=color_normals, showsegments=true)

# Compute Neighbours
# m = computeNeighbours(points, faces)

# is_colored = fill(0, length(faces))
# Red = RGB(1.0,0.0,0.0)
# Blue = RGB(0.0,0.0,1.0)
# colors_neigh = fill(Red, length(faces))

# colors_neigh[1] = Blue # Face 1 goes Blue 
# is_colored[1] = 1

# # 2-Color the Mesh (DFS of connectivity graph)
# for i = 2:length(faces)
#     if is_colored[i] == 1
#         continue
#     end
#     for v in m[i]
#         if colors_neigh[i] == Blue 
#             colors_neigh[v] = Red
#             is_colored[v] = 1
#         elseif colors_neigh[i] == Red
#             colors_neigh[v] = Blue
#             is_colored[v] = 1
#         end
#     end
# end
# viz(mesh, showsegments=true, color=colors_neigh)