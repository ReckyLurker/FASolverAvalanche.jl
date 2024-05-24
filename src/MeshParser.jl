using Meshes, ProgressBars


# Function to Parse Points [allows Multiple Threads]
function parsePoints(location::String)
    t_start = time()
    
    f_pointCloud = open(location, "r")
    s_pointsCloud = split(read(f_pointCloud, String), '\n')
    close(f_pointCloud)
    
    # [Advice on Stack] Only use Regex if you have to
    nPoints = parse(Int, s_pointsCloud[20])
    s_pointsCloud = s_pointsCloud[22:22+nPoints-1]

    println("Reading Point Cloud...")
    pbar_points = ProgressBar(total=length(s_pointsCloud))
    Points = Vector{Point3}(undef, length(s_pointsCloud)) # Preallocating the array
    Threads.@threads for i in eachindex(s_pointsCloud)
        s = split(strip(s_pointsCloud[i], ('(', ')')))
        Points[i] = Point3(parse(Float64, s[1]),parse(Float64, s[2]), parse(Float64, s[3]))
        update(pbar_points)
    end
    println("Parsing Points took: ", time() - t_start, " seconds")
    return Points
end

# Function to parse Faces [Allows Multiple Threads]
function parseFaces(location)
    t_start = time()
    f_faces = open(location, "r")
    s_faces = read(f_faces, String)
    close(f_faces)

    r_nfaces = r"\n\d+"
    nFaces = parse(Int, match(r_nfaces, s_faces).match)

    r_faces = r"\n\d+\(((?:\d+\s*)+)\)"
    s_faces = [s.match for s in eachmatch(r_faces, s_faces)]
    

    faces = Vector(undef, length(s_faces))
    println("Reading Face Information...")
    pbar_faces = ProgressBar(total=length(s_faces))
    Threads.@threads for i in eachindex(s_faces)
        r_pattern = r"\d+"
        res_pattern = eachmatch(r_pattern, s_faces[i])
        pointValues = [parse(Int, p.match)+1 for p in res_pattern]
        faces[i] = tuple(pointValues[2:end] ...)
        update(pbar_faces)
    end
    println("Parsing Faces took: ", time() - t_start, " seconds")
    return faces, nFaces
end

# Function to Parse Face Labels
function parseFaceLabels(locationFaceLabels, Faces, correctionFactor)
    t_start = time()
    f_fafaces = open(locationFaceLabels, "r")
    s_fafaces = read(f_fafaces, String)
    close(f_fafaces)

    s_fafaces = split(s_fafaces, '\n')
    nfaFaces = parse(Int, s_fafaces[20])
    s_fafaces = s_fafaces[22:22+nfaFaces-1]
    println("Reading Face Labels...")
    i_fafaces = [parse(Int, p)+1+correctionFactor for p in s_fafaces]
    fafaces = Faces[i_fafaces]
    println("Parsing Face Labels took: ", time() - t_start, " seconds")
    return fafaces
end


