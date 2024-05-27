include("./MeshParser.jl") # Parses Mesh
include("Mesher.jl") # Pre-Processing Algorithms and cell definitions
import DelimitedFiles.writedlm as write

function generateMesh(locationPoints::String, locationFaces::String, locationFaceLabels::String)
    println("Using ", Threads.nthreads(), " threads")
    Points = parsePoints(locationPoints)
    Faces, nFaces = parseFaces(locationFaces)
    FaFaces = parseFaceLabels(locationFaceLabels,Faces, length(Faces) - nFaces)
    # CorrectionFactor = dev fix, don't change!

    # Remove Unneccessary Points
    idx_FaPoints = [x[j] for x in FaFaces for j in eachindex(x)] # Get the indcies of all points used for FA Scheme
    unique!(idx_FaPoints) # Eliminate Repeated
    sort!(idx_FaPoints) # Sort in case of unordered
    idx_map = Dict(idx_FaPoints .=> 1:length(idx_FaPoints)) # ReMap 
    FaFaces_new = Vector(undef, length(FaFaces))
    FaPoints = Points[idx_FaPoints]
    Threads.@threads for i in eachindex(FaFaces)
        v = [idx_map[j] for j in FaFaces[i]]
        FaFaces_new[i] = tuple(v...)
    end
    return FaPoints, FaFaces_new
end

function generateMesh(points::Vector, faces::Vector)
    # Generate Mesh for visualization    
    connec = connect.(faces, Ngon)
    mesh = SimpleMesh(points, connec)
    return mesh
end


# Cache Optimization needed! - Store extensive computational data [Neighbours, etc.] and reload them if no changes are found in the input file.
function preProcess(locationPoints::String, locationFaces::String, locationFaceLabels::String)
    points, faces = generateMesh(locationPoints, locationFaces, locationFaceLabels)
    normals = calculateNormals(points, faces)
    centers = [p[1] for p in normals] # Face Centers
    edgeCenters = computeEdgeCenters(points, faces) # Compute Edge Centers 
    directions_face = [p[2] for p in normals] # Face Normal directions
    if isfile("./neighbours.txt")
        neighbours = parseNeighbours("./neighbours.txt")
    else 
        neighbours = computeNeighbours(points, faces) # Neighbours of each cell [Cache This by writing to disk]
        
    end
    areas = computeFaceAreas(points, faces, centers) # Areas of each face
    edgeLengths = computeEdgeLengths(points, faces) # Edge Lengths of each edge of each face
    dets = calculateDeterminants(points, faces) # Integral (Jacobian) / Coordinate (DCM) transform determinant
    # Construct Mesh Cells [with computed Data]
    meshCells = Vector(undef, length(faces))
    Threads.@threads for j in eachindex(faces)
        meshCells[j] = Cell(centers[j], j, [faces[j]...], edgeCenters[j], directions_face[j], areas[j], edgeLengths[j],neighbours[j], dets[j], Vec3(0.0,0.0,0.0), 0.0, 0.0)
    end
    return meshCells, points, faces
end

function MeshBounds(Cells)
    x_min = Cells[1].center_coords.coords[1]
    y_min = Cells[1].center_coords.coords[2]
    x_max = Cells[1].center_coords.coords[1]
    y_max = Cells[1].center_coords.coords[2]
    for i in eachindex(Cells)
        coords = Cells[i].center_coords
        x = coords.coords[1]
        y = coords.coords[2]
        x_min = min(x_min, x)
        x_max = max(x_max, x)
        y_min = min(y_min, y)
        y_max = max(y_max, y)
    end
    return [x_min, x_max, y_min, y_max]
end
