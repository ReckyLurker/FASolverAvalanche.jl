import FASolverAvalanche as solver 
import GLMakie as Mke # Visualizer
import Images.RGB as RGB # Color Schemes
import Images.Gray as Gray 
using UnicodePlots
import BenchmarkTools

# Necessary! --- These should be defined in global domain [preferrably non-mutable]
const alpha = 0.5
const zeta = 1.25
const rho = 1.4

# Re-structured Code 

Cells, surface_points, surface_faces = solver.preProcess("/home/tj/Desktop/FASolverAvalanche.jl/data/wolfsgrube_mixed/points",
                    "/home/tj/Desktop/FASolverAvalanche.jl/data/wolfsgrube_mixed/faces",
                    "/home/tj/Desktop/FASolverAvalanche.jl/data/wolfsgrube_mixed/faceLabels")

solver.meshBounds(Cells) # of cell centers 
solver.meshBounds(surface_points) # of vertices 


# # Initialize geometry 
# solver.initializeGeometry([1000.0, 1300.0, 314.0,600.0], Cells, h0 = 0.1) # for wolfsgrube
# solver.initializeGeometry([1700.0, 2200.0, 100.0,700.0], Cells, h0 = 0.1) # for wolfsgrube_mixed

# Initialize Pressure 
solver.initializePressure!(Cells, alpha, zeta, rho, points)

# # Visualize 
# Red = RGB(1.0,0.0,0.0)
# Yellow = RGB(1.0,1.0,0.0)
# Blue = RGB(0.0,0.0,1.0)

# colors = fill(Blue, length(Cells))

# for idx in eachindex(Cells)
#     if Cells[idx].h == 0.1
#         colors[idx] = Red 
#     elseif Cells[idx].pb != 0
#         colors[idx] = Yellow 
#     end
# end

# mesh = solver.createVizMesh(surface_points, surface_faces)
# Mke.plot(mesh, showsegments = true, color=colors)







# OLD CODE [Use OLD Files to use]

# setInitialConditionsPolygon([1000.0, 1300.0, 314.0,600.0], Cells, h0 = 0.1)

# function initpressureg(x)
#     p0 = InitPressure(Cells, fdiff.value(x[1]), fdiff.value(x[2]), fdiff.value(x[3]), points, atol=1e-2)
#     return p0
# end

# initpressureg([0.5, 1.25, 1.4])


# h0 = 0.1 -> Red, else Blue 
# colors = Vector(undef, length(faces))
# Red = RGB(1.0,0.0,0.0)
# Blue = RGB(0.0,0.0,1.0)
# Yellow = RGB(1.0,1.0,0.0)
# # Visualizer
# colors = fill(Blue, length(Cells))
# Threads.@threads for i in eachindex(Cells)
#     if Cells[i].h == 0.1
#         colors[i] = Red
#     elseif p0[i] != 0
#         colors[i] = Yellow 
#     else 
#         colors[i] = Blue
#     end
# end

# mesh = generateMesh(points, faces)
# Mke.plot(mesh, color=colors, showsegments=true, label="Mesh")

# Mesh Generation and visualization


# # Linear Map from range(a,b) to range(c,d) [Assuming x lies in (a,b)]
# function map(a,b,c,d,x)
#     return c + (d - c)/(b-a)*(x-a)
# end

# points, faces = generateMesh("/home/recklurker/RWTHIntern/Julia/points",
#                     "/home/recklurker/RWTHIntern/Julia/faces",
#                     "/home/recklurker/RWTHIntern/Julia/faceLabels")
# mesh = generateMesh(points, faces)
# # To Visualize mesh
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

