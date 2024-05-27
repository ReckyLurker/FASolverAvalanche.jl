import Meshes as Mesher
import LinearAlgebra as LinAlg  
import ProgressBars as bar
import Base.+, Base.*, Base.one

function divide(Center::Vector{Float64}, n::Int)::Mesher.Point3
    return Mesher.Point3(Center[1]/n, Center[2]/n, Center[3]/n)
end

function calculateUnitNormal(vectorA::Mesher.Vec3, vectorB::Mesher.Vec3)::Mesher.Vec3
    return LinAlg.normalize(LinAlg.cross(vectorA, vectorB))
end

function calculateCenter(vertices::Vector{Mesher.Point3})::Mesher.Point3
    center = [0,0,0]
    for vertex in vertices
        center = [sum(x) for x in zip(center, vertex.coords)]
    end
    return divide(center,length(vertices))
end

# Calculate Normals of a Mesh about the center [Allows Multiple Threads]
function calculateNormals(points, faces)::Vector{Tuple{Point3, Vec3}}
    t_start = time()
    println("Using ", Threads.nthreads(), " threads")
    result = Vector(undef, length(faces))
    pbar_normals = ProgressBar(total=length(faces))
    println("Calculating Normals...")
    Threads.@threads for i in eachindex(faces)
        i_vertices = [p for p in faces[i]] # Convert to vector
        vertices = points[i_vertices]
        center = calculateCenter(vertices)
        unitNormal = calculateUnitNormal(vertices[2] - vertices[1], vertices[3] - vertices[1])    
        result[i] = (center, unitNormal)
        update(pbar_normals)
    end
    println("Calculating Normals took: ", time() - t_start, " seconds")
    return result
end

# Calculate Global to Local Transformation determinants
function calculateDeterminants(points, faces)::Vector{Float64}
    t_start = time()
    println("Using ", Threads.nthreads(), " threads")
    result = Vector(undef, length(faces))
    pbar_det = ProgressBar(total=length(faces))
    println("Calculating determinants...")
    Threads.@threads for i in eachindex(faces)
        i_vertices = [p for p in faces[i]] # Convert to vector from tuple
        vertices = points[i_vertices]
        unitNormal = calculateUnitNormal(vertices[2] - vertices[1], vertices[3] - vertices[1])    
        vA = LinAlg.normalize(vertices[2] - vertices[1])
        vB = -calculateUnitNormal(vA, unitNormal)
        result[i] = LinAlg.det(computeTransformationMatrix([vA, vB, unitNormal]))
        update(pbar_det)
    end
    println("Calculating determinants took: ", time() - t_start, " seconds")
    return result
end

# assuming a Global Array of Points and Faces is maintained
mutable struct Cell
    center_coords::Mesher.Point3 # Global Coordinates of the Face Center
    center_idx::Int # Index of the Face
    vertices_idx::Vector{Int} # Indicies of the vertices [For accessing from global points array]
    edgeCenters::Vector{Mesher.Point3}
    faceNormal::Mesher.Vec3 # Surface Normal
    area::Float64 # Area of the Cell
    edgeLengths::Vector{Float64} # Edge Lengths of the edges [same order as of edge list]
    neighbours::Vector{Int} # Index of neighbouring face centers of each point [same order as of edge list]
    det::Float64 # Integral/Coordinate Transform determinant
    vel::Mesher.Vec3 # velocity of fluid [Global Coordinate System]
    h::Float64 # Flow thickness on this Cell
    p::Float64 # Basal Pressure on this Cell
    # Det(J) - Store the determinant of Jacobian (Coordinate Transform), will have to check if its faster to store than to compute at run time
    # localCoords - Array of 3Vectors defining the local coordinate system [will have to check if its faster to store than to compute at run time]
end

# Function to go from LocalCoords to GlobalCoords
function computeTransformationMatrix(localCoords)::Matrix{Float64}
    DCM = Matrix{Float64}(undef, 3, 3)
    globalCoords = eltype(DCM).(one(DCM))
    for i in 1:3
        for j in 1:3
            DCM[i,j] = LinAlg.dot(localCoords[i], globalCoords[:,j])
        end
    end
    return LinAlg.transpose(DCM)
end

# Current Complexity O(pf+e) - Optimization required!
function computeNeighbours(points, faces)
    t_start = time()
    println("Using ", Threads.nthreads(), " threads")
    point_face_map = [Vector{Int}() for _ in eachindex(points)]
    pbar_neigh = ProgressBar(total=length(points))
    println("Computing PointFaceMap...")
    Threads.@threads for k in 1:(length(points)*length(faces)-1)
        i = k ÷ length(faces) + 1
        j = k % length(faces) + 1
        for v in faces[j]
            if i == v 
                push!(point_face_map[i], j)
            end
        end
        if j == 1
            update(pbar_neigh)
        end
    end
    t1 = time()
    println("Computing PointFaceMap took: ", t1 - t_start, " seconds")
    println("Computing Neighbours...")
    neighbours = [Set{Int}() for _ in eachindex(faces)]
    Threads.@threads for j in eachindex(faces)
        for v in eachindex(faces[j])
            e = (v) % length(faces[j]) + 1
            sv = Set(point_face_map[faces[j][v]])
            se = Set(point_face_map[faces[j][e]])
            neig = intersect(sv,se)
            union!(neighbours[j], neig)
            setdiff!(neighbours[j], Set([j]))
        end
    end
    neighbours_corr = [Vector{Int}() for _ in eachindex(faces)]
    for j in eachindex(neighbours)
        neighbours_corr[j] = [neighbours[j]...]
    end
    println("Computing Neighbours took: ", time() - t1, " seconds")

    # Save the computed data to a file!
    neighbours_corr
end

function computeFaceAreas(points, faces, centers)
    areas = Vector(undef, length(faces))
    println("Computing Face Areas...")
    Threads.@threads for i in eachindex(faces)
        area = 0.0
        center = centers[i]
        for j in eachindex(faces[i])
            e = (j) % length(faces[i]) + 1
            area = area + LinAlg.norm2(LinAlg.cross(center - points[faces[i][j]], center - points[faces[i][e]]))*0.5
        end
        areas[i] = area
    end
    return areas
end

function computeEdgeLengths(points, faces)
    edgeLengths = [Vector{Float64}() for _ in eachindex(faces)]
    println("Computing Edge Lengths...")
    Threads.@threads for i in eachindex(faces)
        for j in eachindex(faces[i])
            e = (j) % length(faces[i]) + 1 
            push!(edgeLengths[i], LinAlg.norm2(points[faces[i][e]] - points[faces[i][j]]))            
        end
    end
    return edgeLengths
end

function computeEdgeCenters(points, faces)
    edgeCenters = [Vector{Point3}() for _ in eachindex(faces)]
    Threads.@threads for j in eachindex(faces) 
        for i in eachindex(faces[j])
            ii = i % length(faces[j]) + 1
            push!(edgeCenters[j], 0.5 * (points[faces[j][i]] + points[faces[j][ii]]))
        end
    end
    return edgeCenters
end

# Common Operations
function +(A::Mesher.Point3, B::Mesher.Point3)::Mesher.Point3
    return Mesher.Point3(A.coords .+ B.coords)
end

function *(A::Float64, B::Mesher.Point3)::Mesher.Point3
    return Mesher.Point3(A .* B.coords)
end

function *(A::Mesher.Point3, B::Float64)::Mesher.Point3
    return Mesher.Point3(B .* A.coords)
end

