using QuadGK

"""
Functions required to compute the RRR counts will go in here.
Formulas are from https://mathworld.wolfram.com/Circle-CircleIntersection.html
"""

function A(R, d)
    R^2*acos(d/R) - d*sqrt(R^2 - d^2)
end

"""
    intersect area of two spheres with radii R1 and R2 distance d away from each other.
"""
function intersect(R1, R2, d)
    d1 = (d^2 + R1^2 - R2^2)/2/d
    d2 = (d^2 - R1^2 + R2^2)/2/d
    A(R1, d1) + A(R2, d2)
end

"""
    Volume of intersect of two bins. r3-fixed, r1/r2 vary within a bin.
"""
function fixedR1(r1min, r1max, r2min, r2max, r3)
    intersect(r1max,r2max,r3) - intersect(r1max,r2min,r3) - intersect(r1min,r2max,r3) + intersect(r1min,r2min,r3)
end

"""
    Volume but now integrated over r1 bin.
"""
function RRRVol(r1min, r1max, r2min, r2max, r3min, r3max)
    quadgk(x -> fixedR1(r1min, r1max, r2min, r2max, x)*2*π*x, r3min, r3max, rtol=1e-3)[1]
end

"""
    function RRR(rbinedges)

    Compute the volume available to a triplet (for wp3 computation).
    This will need to be multiplied by N^3/V^2, where N is the number of particles and V the cube volume.

    Input:
    rbinedges - array of edges of the bins.
"""
function RRR(rbinedges)
    rmid = (rbinedges[1:end-1] + rbinedges[2:end])/2
    dr = rbinedges[2] - rbinedges[1]
    Nbins = length(rmid)

    Ncol = 0
    for i1 in 1:Nbins, i2 in i1:Nbins, i3 in i2:Nbins
        r1 = rmid[i1] # The shortest
        r2 = rmid[i2]
        r3 = rmid[i3] # The longest
        if r3 + dr/2 < r2 - dr/2 + r1 - dr/2
            Ncol += 1
        end
    end

    hist_RRR = zeros(Ncol, 4)
    Ncol = 1
    for i1 in 1:Nbins, i2 in i1:Nbins, i3 in i2:Nbins
        r1 = rmid[i1] # The shortest
        r2 = rmid[i2]
        r3 = rmid[i3] # The longest
        if r3 + dr/2 < r2 - dr/2 + r1 - dr/2
            hist_RRR[Ncol, 1:3] = [r1 r2 r3]
            hist_RRR[Ncol, 4] = RRRVol(rbinedges[i1], rbinedges[i1+1], rbinedges[i2], rbinedges[i2+1], rbinedges[i3], rbinedges[i3+1])
            if r1 == r2 == r3
                2 + 2
            elseif r1 == r2 || r2 == r3
                hist_RRR[Ncol, 4] *= 3
            else
                hist_RRR[Ncol, 4] *= 6
            end
            Ncol += 1
        end
    end
    return hist_RRR
end
