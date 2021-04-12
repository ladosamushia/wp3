include("wp3.jl")

"""
    neighbouring_pairs()

    Output:
    - aa: An array of quadrupoles of integers.

    In 2D you have 9 neighbours (including yourself). 
"""
function neighbouring_pairs()
    aa = Array{Int64,1}[]
    for i in -1:1, j in -1:1
        push!(aa, [i, j])
    end
    return aa
end

"""
    histogram_dd!(xy1, xy2, dr, hist)

    look at all pairs in 2D arrays xy1, xy2.

    Input:
    xy? - Float 2D array.
    dr - Float. Bin width.
    hist - Int 3D array.

    xy? can be identical, in which case this computes auto-triplets, but does not divide by the simmetry factor.
"""

function histogram_dd!(xy1, xy2, dr, hist)

    N = size(hist)[1]
    # Avoid double counting identical cells
    for i1 in 1:size(xy1)[2], i2 in 1:size(xy2)[2]
        p1 = view(xy1, :, i1)
        p2 = view(xy2, :, i2)
        r12 = distance(p1, p2)
        if r12 == 0 continue end
        h = ceil(Int, r12/dr)
        if h > N 
            continue
        else
            hist[h] += 1
        end
    end
    
end

"""
    double_loop(xy_cube1, xy_cube2, Ngal1, Ngal2, dr, hist)

    Loop over all pairs of neighbours in the xy_cube grid and histogram them.

    Input:
    xy_cube? - 4D float array of x y coordinates. The first two indeces are for the x and y, the last two are indexing the cell.
    Ngal - 2D Int array of how many galaxies ended up in each cell.
    dr - Float. Bin width.
    hist - 1D Int array of histogrammed separations.

    Looks at all pairs that are separatd by no more than one cell and histograms them.
"""
function double_loop!(xy_cube1, xy_cube2, Ngal1, Ngal2, dr, hist)

    Ncube = size(xy_cube1)[end] - 2
    ss = neighbouring_pairs()
    for ix in 2:Ncube+1, iy in 2:Ncube+1
        for (jx, jy) in ss
        xy1 = view(xy_cube1, :, 1:Ngal1[ix, iy], ix, iy)
        xy2 = view(xy_cube2, :, 1:Ngal2[ix+jx, iy+jy], ix+jx, iy+jy)
        histogram_dd!(xy1, xy2, dr, hist)
        end
    end
end

"""
    DDD(x, y, dr, Nbin)

    Input:
    x - Array of floats
    y - Array of floats
    dr - Float. wp3 bin size.
    Nbin - Int. Number of bins in each r.

    Output:
    rhist - 2D array. Floats. The columns are r1, r2, r3, wp3.
"""
function DD(x1, y1, x2, y2, dr, Nbin)
    L = maximum([x1 x2 y1 y2]) - minimum([x1 x2 y1 y2])
    Ncell = floor(Int, L/(dr*Nbin))
    xy1_cube, N1_cube = make_cube(x1, y1, Ncell)
    xy2_cube, N2_cube = make_cube(x2, y2, Ncell)
    hist = zeros(Nbin)
    double_loop!(xy1_cube, xy2_cube, N1_cube, N2_cube, dr, hist)
    # Account for self-paris
    if x1 == x2
        hist ./= 2
    end
    return hist
end