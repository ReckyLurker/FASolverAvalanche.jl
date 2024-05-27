import LinearAlgebra.norm2 as magnitude
import LinearAlgebra as LinAlg
import Meshes as Mesher
import ForwardDiff as fdiff 

const g = Mesher.Vec3(0.0,0.0,-9.81)

# For scalar second order central interpolation [velocity in global coords, scalar variables]
@inline function interpolate(center_face1, center_face2, edgeCenter, vars)
    pe = edgeCenter - center_face1
    en = center_face2 - edgeCenter
    pen = pe + en
    mpe = magnitude(pe)
    mpen = magnitude(pen)
    frac = mpe/mpen
    return frac*vars[1] + (1 - frac)*vars[2]
    # 0.25 * (vars[1])  + 0.75*(vars[2])
end


# Calculate Hessians for Edge Contributions O(E*T(hessian))
function computeHessians(Cells; func=interpolate)
    hessians = [Vector() for _ in eachindex(Cells)]
    for i in eachindex(Cells)
        CellP = Cells[i]
        j = 0
        for idx in CellP.neighbours
            CellE = Cells[idx]
            j = j + 1
            function pContr(x)
                l = div(length(x), 2)
                hvalues = x[1:l]
                pvalues = x[l+1:end]
                h = func(CellP.center_coords, CellE.center_coords, CellP.edgeCenters[j], hvalues)
                p = func(CellP.center_coords, CellE.center_coords, CellP.edgeCenters[j], pvalues)
                return h*p
            end
            hessian_edge = fdiff.hessian(pContr, [CellP.h, CellE.h, CellP.p, CellE.p])
            push!(hessians[i], hessian_edge)
        end
    end
    return hessians
end


# Test Membership of a test point and a bounding region
function testInside(coords, boundingCoords)::Bool
    angle_made = 0.0
    for i in eachindex(boundingCoords)
        e = i % length(boundingCoords) + 1
        vA = Point2(coords.coords[1:2]) - boundingCoords[i]
        vB = Point2(coords.coords[1:2]) - boundingCoords[e]
        lvA = magnitude(vA)
        lvB = magnitude(vB)
        angle_made = angle_made + acos(LinAlg.dot(vA,vB)/(lvA*lvB))
    end
    return abs(angle_made - 2*pi) < 0.1 ? true : false
end

# The Points given bound a polygon: coord_bounds contain the points -> [x_min, x_max, y_min, y_max]
function setInitialConditionsPolygon(coord_bounds, Cells; h0 = 0.0, u0=Mesher.Vec3(0.0,0.0,0.0), epsilon = 1e-6)
    insideCells = BoundingPolygon(Cells, coord_bounds, epsilon = epsilon)
    # Set Initial Conditions
    for cell in insideCells
        cell.h = h0
        cell.vel = u0
    end
end

# Find all the points inside the region in the range x (xmin, xmax) and y (ymin, ymax)
# Choose Epsilon as required accuracy, default = 1e-6
function BoundingPolygon(Cells, coord_bounds; epsilon = 1e-6)
    x_min = coord_bounds[1] - epsilon
    x_max = coord_bounds[2] + epsilon
    y_min = coord_bounds[3] - epsilon
    y_max = coord_bounds[4] + epsilon
    boundingCoords = [Point2(x_min, y_min), Point2(x_min, y_max), Point2(x_max, y_max), Point2(x_max, y_min)]
    insideCells = []
    for cell in Cells
        if(testInside(cell.center_coords, boundingCoords))
            push!(insideCells, cell)
        end
    end
    return insideCells
end

# Need Symbolics for General Solvers  
function InitPressure(Cells, alpha, zeta, rho)
    A = zeros(length(Cells), length(Cells))
    B = zeros(length(Cells))
    rho_i = 1/rho 
    hessians = computeHessians(Cells)
    # For each cell 
    Threads.@threads for i in eachindex(Cells)
        A[i,i] = alpha * rho_i * Cells[i].area
        B[i] = Cells[i].det * LinAlg.dot(Cells[i].faceNormal, g) * Cells[i].h * Cells[i].area
        for j in eachindex(Cells[i].neighbours)
            n = Cells[i].neighbours[j]
            A[i,i] += (hessians[i][j][3,1] * Cells[i].h + hessians[i][j][4,2] * Cells[n].h) * alpha * Cells[i].det * rho_i * LinAlg.dot(Cells[i].faceNormal, Cells[n].center_coords - Cells[i].center_coords) * Cells[i].edgeLengths[j]
            A[i,n] = alpha * Cells[i].det * rho_i * LinAlg.dot(Cells[i].faceNormal, Cells[n].center_coords - Cells[i].center_coords) * Cells[i].edgeLengths[j] * (hessians[i][j][4,1] * Cells[n].h + hessians[i][j][3,2] * Cells[i].h)
            
            function uContr(hvalues, uvalues)
                h = interpolate(Cells[i].center_coords, Cells[n].center_coords, Cells[i].edgeCenters[j], hvalues)
                u = interpolate(Cells[i].center_coords, Cells[n].center_coords, Cells[i].edgeCenters[j], uvalues)
                return zeta * Cells[i].det * LinAlg.dot(Cells[i].faceNormal, u) * LinAlg.dot(Cells[n].center_coords - Cells[i].center_coords, u) * h * Cells[i].edgeLengths[j]
            end
            B[i] += uContr([Cells[i].h, Cells[n].h], [Cells[i].vel, Cells[n].vel])
        end
    end
    A\B 
end




