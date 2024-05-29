include("./meshParser.jl") # For Parsing input files 
import JLD2, Meshes 

# Parse mesh into memory from generated text files 
function parseMesh(location_points::String, location_faces::String, location_face_labels::String)::Tuple{Vector{Meshes.Point}, Vector{Vector{Int}}}
    println("Using ", Threads.nthreads(), " threads")
    mesh_points = parsePoints(location_points)
    mesh_faces, total_faces = parseFaces(location_faces)
    surface_aligned_faces = parseFaceLabels(location_face_labels, mesh_faces, length(mesh_faces) - total_faces)
    # Correction Factor, dev fix, don't change!
    
    # Remove Unneccessary
    surface_points_idx = [x[j] for x in surface_aligned_faces for j in eachindex(x)] # Get points on the surface aligned mesh 
    unique!(surface_points_idx) # Remove repeated indcies
    sort!(surface_points_idx) # Sort the indices, prepare for remapping 
    
    # Remapping 
    global_surface_idx_map = Dict(surface_points_idx .=> 1:length(surface_points_idx)) # New Index ReMap
    surface_aligned_faces_remapped = Vector(undef, length(surface_aligned_faces))
    surface_points = mesh_points[surface_points_idx]
    Threads.@threads for i in eachindex(surface_aligned_faces)
        surface_aligned_faces_remapped[i] = [global_surface_idx_map[j] for j in surface_aligned_faces[i]]
    end
    return surface_points, surface_aligned_faces_remapped
end

# Generate mesh for visualization
function createVizMesh(mesh_points::Vector{Meshes.Point}, mesh_faces::Vector{Vector{Int}})
    # Tuple the Vector of vectors 
    mesh_faces_tuple = Vector(undef, length(mesh_faces))
    Threads.@threads for i in eachindex(mesh_faces)
        mesh_faces_tuple[i] = tuple(mesh_faces[i]...)
    end

    connec = Meshes.connect.(mesh_faces_tuple, Meshes.Ngon) # Connectivity List 
    mesh = Meshes.SimpleMesh(mesh_points, connec) # Mesh 
end

function preProcess(location_points::String, location_faces::String, location_face_labels::String)
    surface_points, surface_aligned_faces = parseMesh(location_points, location_faces, location_face_labels) # Parse Mesh to Memory 
    surface_aligned_face_normals = computeNormals(surface_points, surface_aligned_faces) # Compute Face Normals 
    surface_aligned_face_centers = [p[1] for p in surface_aligned_face_normals] # Face Centers 
    surface_aligned_face_normal_directions = [p[2] for p in surface_aligned_face_normals] # Normal Directions 
    surface_aligned_edge_centers = computeEdgeCenters(surface_points, surface_aligned_faces) # Edge Centers 

    # Cache Neighbours [Expensive Computation]
    if Base.isfile("../stored/neighbours.jld2")
        surface_aligned_neighbours = JLD2.load_object("../stored/neighbours.jld2")
    else 
        if !Base.isdir("../stored")
            mkdir("../stored")
        else 
            Base.rm("../stored", recursive=true)
            Base.mkdir("../stored")
        end
        surface_aligned_neighbours = computeNeighbours(surface_aligned_points, surface_aligned_faces) # Neighbours for each edge of a face
        JLD2.save_object("../stored/neighbours.jld2")
    end
    surface_aligned_face_areas = computeFaceAreas(surface_points, surface_aligned_points, surface_aligned_face_centers) # Surface Area of each face 
    surface_aligned_face_edge_lengths = computeEdgeLengths(surface_points, surface_aligned_faces) # Edge Lengths of each edge of a face 
    surface_aligned_jacobian_determinant = computeDeterminants(surface_points, surface_faces) # Integral Transform (Jacobian)/ Coordinate Transform Matrix Determinant 

    # Construct Mesh Cells [Initialize to Zero]
    surface_aligned_mesh_cells = Vector(undef, length(surface_aligned_faces))
    Threads.@threads for i in eachindex(surface_aligned_faces)
        surface_aligned_mesh_cells[i] = Cell(surface_aligned_face_centers[i], i, surface_aligned_faces[i], surface_aligned_edge_centers[i], surface_aligned_face_normal_directions[i], surface_aligned_face_areas[i], surface_aligned_face_edge_lengths[i], surface_aligned_neighbours[i], surface_aligned_jacobian_determinant[i], Meshes.Vec3(0.0,0.0,0.0), 0.0, 0.0)
    end
    return surface_aligned_mesh_cells, surface_points, surface_aligned_faces
end