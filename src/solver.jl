import LinearAlgebra.norm2 as magnitude

# For scalar second order central interpolation [velocity in global coords, scalar variables]
function interpolate(center_face1, center_face2, edgeCenter, var_face1, var_face2)
    pe = edgeCenter - center_face1
    en = center_face2 - edgeCenter
    pen = pe + en
    mpe = magnitude(pe)
    mpen = magnitude(pen)
    frac = mpe/mpen
    return frac*var_face1 + (1 - frac)*var_face2
end



# Simulation Parameters
zeta = 1.25 # Velocity Depth Profile Assumption
alpha = 0.5 # Pressure Depth Profile Assumption


