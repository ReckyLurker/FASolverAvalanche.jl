import Meshes.Point, Meshes.Vec 
import LinearAlgebra as LinAlg 
import ProgressBars as bar 
import Base.+, Base.*


############################################
#              CELL STRUCTURE              #
############################################

# For each face 
mutable struct Cell 
    center_coords::Point # Coords of the face center 
    center_idx::Int # Index of the Cell 
    vertices::Vector{Point} # Coords of the vertices
    vertices_idx::Vector{Int} # Indicies of the vertices 
    edge_centers::Vector{Point} # Coords of edge centers [same order as of edges]
    edge_lengths::Vector{Float64} # Lengths of edges [same order as edges]
    face_normal::Vec # Surface Normal 
    area::Float64 # Area of the Cell 
    neighbours::Vector{Int} # Indicies of faces sharing an edge 
    det::Float64 # Integral Transform/Jacobian/Coordinate transform matrix's determinant
    vel::Vec # Velocity of shallow fluid [global coordinates]
    h::Float64 # Snow Thickness on this cell 
    pb::Float64 # Basal Pressure on this cell 
end

############################################
#              UTILITY FUNCTIONS           #
############################################

# Requires closed polygon checks 
# Returns local coordinate system of the cell [3 orthogonal vectors] 
function localCoords(cell::Cell)::Vector{Vec}
    idx1 = rand(1:length(cell.vertices))
    idx2 = idx1 % length(cell.vertices) + 1
    idx3 = idx2 % length(cell.vertices) + 1
    edge1 = LinAlg.normalize(cell.vertices[idx2] - cell.vertices[idx1])
    edge2 = LinAlg.normalize(cell.vertices[idx3] - cell.vertices[idx2])
    surface_normal = normalized_cross(edge1, edge2)
    bi_normal = -normalized_cross(edge1, surface_normal)
    return [edge1, bi_normal, surface_normal]
end

# Requires closed polygon checks
# Returns the local coordinate system for a closed polygonal face.
function localCoords(vertices::Vector{Point})
    l = length(vertices)
    if l < 3
        throw("required 3 vertices, got $l")
    end
    idx1 = rand(1:l) 
    idx2 = idx1 % l + 1
    idx3 = idx2 % l + 1
    edge1 = LinAlg.normalize(vertices[idx2] - vertices[idx1])
    edge2 = LinAlg.normalize(vertices[idx3] - vertices[idx2])
    surface_normal = normalized_cross(edge1, edge2)
    bi_normal = -normalized_cross(edg1, surface_normal)
    return [edge1, bi_normal, surface_normal]
end

# Returns the transformation matrix for transforming local coords to global coords
function computeTransformationMatrix(local_coords::Vector{Vec})::Matrix{Float64}
    DCM = Matrix{Float64}(undef, 3, 3)
    global_coords = eltype(DCM).(one(DCM))
    for i in 1:3
        for j in 1:3
            DCM[i,j] = LinAlg.dot(local_coords[i], global_coords[:,j])
        end
    end
    return LinAlg.transpose(DCM)
end

# Returns normalized cross-product of two vectors
function normalized_cross(vector1::Vec, vector2::Vec)::Vec
    return LinAlg.normalize(LinAlg.cross(vector1, vector2))
end 

# Returns Jacobian/Integral/Coordinate transformation determinants [only for solver_initialization]
function computeDeterminants(points::Vector{Point}, faces::Vector{Vector{Int}})
    t_start = Base.time()
    cell_determinants = Vector(undef, length(faces))
    println("Calculating determinants...")
    Threads.@threads for i in eachindex(faces)
        vertices = points[faces[i]]
        local_coords = localCoords(vertices)
        cell_determinants[i] = LinAlg.det(computeTransformationMatrix(local_coords))
    end
    println("Calculating determinants took: ", Base.time() - t_start, " seconds")
    return cell_determinants
end

# Returns the arithmetic mean of all vertices 
function computeCenter(vertices::Vector{Point})::Point 
    return sum([p.coords for p in vertices]) ./ length(vertices) + Point(0.0,0.0,0.0)
end

# Returns normals of the given geometry [only for solver initialization]
function computeNormals(points::Vector{Point}, faces::Vector{Vector{Int}})::Vector{Tuple{Point, Vec}}
    t_start = Base.time()
    normals_centers_directions = Vector(undef, length(faces))
    println("Calculating normals...")
    Threads.@threads for i in eachindex(faces)
        vertices = points[faces[i]]
        local_coords = localCoords(vertices)
        surface_normal = local_coords[3]
        center = computeCenter(vertices)
        normal_centers_directions[i] = (center, surface_normal)
    end
    println("Computing normals took: ", Base.time() - t_start, " seconds")
    return normals_centers_directions
end

# Needs Optimization O(pf+e) is the current complexity 
# Returns neighbours of each cell for a given collection of cells. [Initialization Only]
function computeNeighbours(points::Vector{Point}, faces::Vector{Vector{Int}})::Vector{Vector{Int}}
    t_start = Base.time()
    point_face_map = [Vector{Int}() for _ in eachindex(points)]
    pbar_neigh = ProgressBars.ProgressBar(total=length(points))
    println("Computing point face map...")
    Threads.@threads for k in 1:(length(points)*length(faces) - 1)
        i = k ÷ length(faces) + 1
        j = k % length(faces) + 1
        if i in faces[j]
            push!(point_face_map[i], j)
        elseif j == 1 
            ProgressBar.update(pbar_neigh)
        end
    end
    t1 = Base.time()
    println("Computing point face map took: ", t1 - t_start, " seconds")
    neighbours = [Set{Int}() for _ in eachindex(faces)]
    Threads.@threads for j in eachindex(faces)
        for v in eachindex(faces[j])
            e = (v) % length(faces[j]) + 1
            set_v = Set(point_face_map[faces[j][v]])
            set_e = Set(point_face_map(faces[j][e]))
            neighbour = intersect(set_v, set_e)
            union!(neighbours[j], neighbour)
            setdiff!(neighbours[j], Set([j]))
        end
    end
    neighbours_arr = [[p...] for p in neighbours]
    println("Computing neighbours took: ", Base.time() - t1, " seconds")
    return neighbours_arr 
end

# Returns the area of a collection of cells [Solver INIT ONLY]
function computeFaceAreas(points::Vector{Point}, faces::Vector{Vector{Int}}, centers::Vector{Point})
    areas = Vector(undef, length(faces))
    println("Computing face areas...")
    Threads.@threads for i in eachindex(faces)
        area = 0.0
        center = centers[i]
        for j in eachindex(faces[i])
            e = j % length(faces[i]) + 1
            internal_edge1 = center - points[faces[i][j]]
            internal_edge2 = center - points[faces[i][e]]   
            area = area + LinAlg.normalize(LinAlg.cross(internal_edge1, internal_edge2)) * 0.5
        end
        areas[i] = area 
    end
    return areas 
end

# Returns edge lengths of a given geometry [Solver INIT ONLY]
function computeEdgeLengths(points::Vector{Point}, faces::Vector{Vector{Int}})::Vector{Vector{Float64}}
    edge_lengths = [Vector{Float64}() for _ in eachindex(faces)]
    println("Computing Edge Lengths...")
    Threads.@threads for i in eachindex(faces)
        for j in eachindex(faces[i])
            e = j % length(faces[i]) + 1
            push!(edge_lengths[i], LinAlg.normalize(points[faces[i][e]] - points[faces[i][j]]))
        end
    end
    return edge_lengths
end

# Returns edge centers for each edge of the given geometry [SOLVER INIT ONLY]
function computeEdgeCenters(points::Vector{Point}, faces::Vector{Vector{Int}})::Vector{Vector{Point}}
    edge_centers = [Vector{Point}() for _ in eachindex(faces)]
    Threads.@threads for j in eachindex(faces)
        for i in eachindex(faces[j])
            ii = i % length(faces[j]) + 1
            push!(edgeCenters[j], 0.5 * (points[faces[j][i]].coords + points[faces[j][ii]].coords + Point(0.0,0.0,0.0)))
        end
    end
    return edge_centers
end


############################################
#              COMMON OPERATORS            #
############################################

function +(A::Point, B::Point)::Point 
    return A.coords .+ B.coords + Point(0.0,0.0,0.0)
end

function *(A::Float64, B::Point)::Point 
    return Point(0.0,0.0,0.0) + A .* B.coords
end

function *(A::Point, B::Float64)::Point 
    return Point(0.0,0.0,0.0) + B .* A.coords 
end
