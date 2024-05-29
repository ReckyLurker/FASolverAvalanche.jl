import Meshes, ProgressBars 

# Function to parse points [on multiple threads]
# Array Size Mutations Present!
function parsePoints(location::String)::Vector{Meshes.Point3}
    t_start = Base.time()
    
    # Read file 
    f_point_cloud = Base.open(location, "r")
    s_points_cloud = Base.split(read(f_point_cloud, String), '\n')
    Base.close(f_point_cloud)

    # [Advice on Stack] Only use Regex if you have to! 
    total_points = Base.parse(Int, s_points_cloud[20])
    s_points_cloud = s_points_cloud[22:22+total_points-1]

    # Process the string
    println("Reading point cloud...")
    pbar_points = ProgressBars.ProgressBar(total=length(s_points_cloud))
    mesh_points = Vector{Meshes.Point}(undef, length(s_points_cloud)) # Preallocating for speed 
    Threads.@threads for i in eachindex(s_points_cloud)
        s = split(strip(s_points_cloud[i], ('(', ')')))
        mesh_points[i] = Meshes.Point3.(parse.(Float64, s))
        ProgressBar.update(pbar_points)
    end
    println("Parsing points took: ", Base.time() - t_start, " seconds")
end

# Function to parse mesh faces
function parseFaces(location::String)::Vector{Vector{Int}}
    t_start = Base.time()

    # Read File 
    f_faces = Base.open(location, "r")
    s_faces = Base.read(f_faces, Base.String)
    Base.close(f_faces)

    # Regex processing 
    pattern_total_faces = r"\n\d+"
    total_faces = Base.parse(Int, match(pattern_total_faces, s_faces).match)

    pattern_faces = r"\n\d+\(((?:\d+\s*)+)\)"
    s_faces_arr = [s.match for s in eachmatch(pattern_faces, s_faces)]

    # String Processing 
    mesh_faces = Vector(undef, length(s_faces_arr))
    println("Reading face infomation...")
    pbar_faces = ProgressBars.ProgressBar(total=length(s_faces_arr))
    Threads.@threads for i in eachindex(s_faces_arr)
        pattern_vertex_idx = r"\d+"
        vertex_idxs = [parse(Int, p.match)+1 for p in eachmatch(pattern_vertex_idx, s_faces_arr[i])]
        faces[i] = vertex_idxs[2:end]
        ProgressBar.update(pbar_faces)
    end
    println("Parsing faces took: ", Base.time() - t_start, " seconds")
    return mesh_faces, total_faces 
end

# Function to parse face labels 
function parseFaceLabels(location::String, mesh_faces::Vector{Vector{Int}}, correction_factor::Int)
    t_start = Base.time()

    # Read file
    f_surface_aligned_faces = Base.open(location, "r")
    s_surface_aligned_faces = Base.split(Base.read(f_surface_aligned_faces, String), '\n')
    Base.close(f_surface_aligned_faces)

    total_surface_aligned_faces = Base.parse(Int, s_surface_aligned_faces[20])
    s_surface_aligned_faces_arr = s_surface_aligned_faces[22:22+total_surface_aligned_faces-1]
    
    # Parse face labels
    println("Reading face labels...")
    idx_surface_aligned_faces = [parse(Int, p)+1+correctionFactor for p in s_surface_aligned_faces_arr]
    surface_aligned_faces = mesh_faces[idx_surface_aligned_faces]
    println("Parsing face labels took: ", Base.time() - t_start, " seconds")
    return surface_aligned_faces
end