import LinearAlgebra.norm2 as magnitude 
import LinearAlgebra as LinAlg 
import LinearSolve as LinSolv  
import ExtendableSparse as ESparse
import Meshes
import ForwardDiff as FDiff 

# Test Membership of a point and a bounding polygon
function testInside(coords::Meshes.Point, bounding_coords::Vector{Meshes.Point2})::Bool
    angle_made = 0.0
    for i in eachindex(bounding_coords)
        e = i % length(bounding_coords) + 1
        vector_a = Meshes.Point2(coords.coords[1:2]) - bounding_coords[i]
        vector_b = Meshes.Point2(coords.coords[1:2]) - bounding_coords[e]
        length_vector_a = magnitude(vector_a)
        length_vector_b = magnitude(vector_b)
        angle_made = angle_made + acos(LinAlg.dot(vector_a, vector_b)/(length_vector_a * length_vector_b))
    end
    return abs(angle_made - 2*pi) < 0.1 ? true : false 
end

# Returns a regular polygon with the given constraints: [x_min, x_max, y_min, y_max]
function findRegularPolygon(coord_limits::Vector{Float64}; npoints = 4)::Vector{Meshes.Point2}
    l = length(coord_limits)
    if l != 4 
        throw("expected 4, got $l")
    end
    x_min = coord_limits[1]
    x_max = coord_limits[2]
    y_min = coord_limits[3]
    y_max = coord_limits[4]
    polygon_center = Meshes.Point2(0.5*(x_min+x_max), 0.5*(y_min+y_max))
    seed_angle = (2*pi)/npoints
    polygon_points = Vector(undef, npoints)
    polygon_points[1] = Meshes.Point2(0.5*(x_min+x_max), y_max) # Initialize inital point
    # Compute rest of the points 
    for i in range(2, npoints)
        vector_to_be_rotated = polygon_points[i-1] - polygon_center
        rotation_matrix = [Base.cos(seed_angle) Base.sin(seed_angle); -Base.sin(seed_angle) Base.cos(seed_angle)] 
        polygon_points[i] = polygon_center + Meshes.Point2(rotation_matrix * vector_to_be_rotated)
    end
    return polygon_points
end

# Find All Cells inside a regular polygon given by the values: [x_min, x_max, y_min, y_max], inclusion of corners by epsilon.
function cellsInsideBoundingPolygon(Cells::Vector{Cell}, coord_limits::Vector{Float64}; npoints = 4, epsilon = 1e-6)::Vector{Int}
    point_uncertainty = [-epsilon, epsilon, -epsilon, epsilon]
    polygon_points = findRegularPolygon(coord_limits .+ point_uncertainty, npoints = npoints)
    cells_inside_polygon = []
    for i in eachindex(Cells)
        if(testInside(Cells[i].center_coords, polygon_points))
            push!(cells_inside_polygon, i)
        end
    end
    return cells_inside_polygon
end

# For scalar second order central interpolation [scalar variables: h and pb]
@inline function central_interpolate(face_center1::Meshes.Point3, face_center2::Meshes.Point3, edge_center::Meshes.Point3, vars::Vector)::eltype(vars[1])
    l = length(vars)
    if l != 2
        throw("Expected 2, got $l")
    end
    pe = edge_center - face_center1
    pen = face_center_2 - face_center1
    frac = magnitude(pe)/magnitude(pen)
    return frac*vars[1] + (1-frac)*vars[2]
end

# For vector second order central interpolation [vector variables: surface_velocity]
function central_interpolate(face_center1::Meshes.Point3, face_center2::Meshes.Point3, edge_center::Meshes.Point3, local_coords1::Vector{Meshes.Vec3}, local_coords2::Vector{Meshes.Vec3}, vars::Vector)::eltype(vars[1])
    l = length(vars)
    if l != 2 
        throw("Expected 2, got $l")
    end
    pe = edge_center - face_center1
    pen = face_center2 - face_center1
    frac = magnitude(pe)/magnitude(pen)
    
    # Transformation Matrices 
    Tp = computeTransformationMatrix(local_coords1)
    Tn = computeTransformationMatrix(local_coords2)
    local_coords_edge = frac .* local_coords1 .+ (1-frac) .* local_coords2
    Te = computeTransformationMatrix(local_coords_edge)

    return LinAlg.transpose(Te) * (frac * Tp * vars[1] + (1-frac) * Tn * vars[2])
end

# Currently limited to 2 values for interpolation at edge, fix this in future!
# Allows user-defined interpolation functions [which are linear, incorrect results may be obtained for discontinuous/non-linear interpolators] 
# Compute Hessians for pressure at edge equations [constraint and momentum eqns]
function computeHessians(Cells; interpolator = central_interpolate)
    hessians = [Vector() for _ in eachindex(Cells)]
    Threads.@threads for i in eachindex(Cells)
        for j in eachindex(Cells[i].neighbours)
            idx = Cells[i].neighbours[j]
            function pressureContribution(x)
                l = div(length(x), 2)
                hvalues = x[1:l]
                pvalues = x[l+1:end]
                h_edge = func(Cells[i].center_coords, Cells[idx].center_coords, Cells[i].edge_centers[j], hvalues)
                p_edge = func(Cells[i].center_coords, Cells[idx].center_coords, Cells[i], edge_centers[j], pvalues)
                return h_edge*p_edge
            end
            hessian_edge = fdiff.hessian(pressureContribution, [Cells[i].h, Cells[idx].h, Cells[i].pb, Cells[idx].pb])
            push!(hessians[i], hessian_edge)
        end
    end
    return hessians
end

# Currently limited to 2 values for interpolation at edge, fix this in future!
# v: current cell, e: neighbour, j: jth edge of current cell, zeta: multiplicative factor
function velocity_contribution(Cells, v, e, j, zeta)
    h_edge = central_interpolate(Cells[v].center_coords, Cells[e].center_coords, Cells[v].edge_centers[j], [Cells[v].h, Cells[e].h])
    u_edge = central_interpolate(Cells[v].center_coords, Cells[e].center_coords, Cells[v].edge_centers[j], localCoords(Cells[v]), localCoords(Cells[e]), [Cells[v].vel, Cells[e].vel])
    return zeta * Cells[v].det * LinAlg.dot(Cells[e].face_normal, u_edge) * LinAlg.dot(Cells[e].center_coords - Cells[v].center_coords, u_edge) * h_edge * Cells[v].edge_lengths[j]
end

# Currently limited to 2 values for interpolation at edge, fix this in future!
# Need symbolics for discontinuos/non-linear interpolators. Using ForwardDiff for general linear interpolators, with user defined solver. 
function initializePressure(Cells, alpha, zeta, rho; atol=1e-5, solver=linsolv.KrylovJL_GMRES())
   A = ESparse.ExtendableSparseMatrix(length(Cells), length(Cells))
   B = zeros(length(Cells))
   rho_inv = 1/rho
   hessians = computeHessians(Cells)

#    For each Cell 
    for i in eachindex(Cells)
        A[i,i] = rho_i * Cells[i].area
        B[i] = Cells[i].det * LinAlg.dot(Cells[i].faceNormal, g) * Cells[i].h * Cells[i].area 
        for j in eachindex(Cells[i].neighbours)
            n = Cells[i].neighbours[j]
            A[i,i] += (hessians[i][j][3,1] * Cells[i].h + hessians[i][j][4,2] * Cells[n].h) * alpha * Cells[i].det * rho_inv * LinAlg.dot(Cells[i].face_normal, Cells[n].center_coords - Cells[i].center_coords) * Cells[i].edge_lengths[j]
            A[i,n] = alpha * Cells[i].det * rho_inv * LinAlg.dot(Cells[i].face_normal, Cells[n].center_coords - Cells[i].center_coords) * Cells[i].edge_lengths[j] * (hessians[i][j][4,1] * Cells[n].h + hessians[i][j][3,2] * Cells[i].h)
            B[i] += velocity_contribution(Cells, i, n, j, zeta)
        end
    end

    # Linear System of Equations for initializing pressure 
    prob = linsolv.LinearProblem(A, B)
    sol = linsolv.solve(prob, solver, abstol=atol, progress=true)
    return sol.u
end