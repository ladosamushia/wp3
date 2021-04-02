"""
    distance(p1, p2)

    Cartesian distance.

    Input:
    p1, p2 - Float{2} contains x and y components.

    Output:
    dist - Float.
"""
function distance(p1, p2)

    dx = (p1[1] - p2[1])
    #dx > L/2 ? dx -= L : dx = dx
    dy = (p1[2] - p2[2])
    #dy > L/2 ? dy -= L : dy = dy
    dist = sqrt(dx.^2 + dy.^2)

end

"""
    neighbouring_triplets()

    Output:
    - aa: An array of quadrupoles of integers.

    In 2D you have 9 neighbours (including yourself). All pairs of neighbours that are also neighbours.
    Returns 4 numbers between -1 and 1 (x and y steps).
"""
function neighbouring_triplets()
    aa = Array{Int64,1}[]
    for i in -1:1, j in -1:1, k in -1:1, l in -1:1
        if abs(i-k) < 2 && abs(j-l) < 2
            push!(aa, [i, j, k, l])
        end
    end
    return aa
end

"""
    hist_index(r1, r2, r3, dr)

    Which bin does this r1, r2, r3 combination belongs to if the bin width is dr?

    Inputs:
    all Floats.

    Output:
    triplet of integers.

    Assumes binning is uniform and starts at 0. Linear binning.
"""
function hist_index(r1, r2, r3, dr)

    ind1 = ceil(Int, r1/dr)
    ind2 = ceil(Int, r2/dr)
    ind3 = ceil(Int, r3/dr)
    return ind1, ind2, ind3

end

"""
    histogram!(xy1, xy2, xy3, dr, hist)

    look at all triplets in 2D arrays xy1, xy2, xy3.

    Input:
    xy? - Float 2D array.
    dr - Float. Bin width.
    hist - Int 3D array.

    xy? can be identical, in which case this computes auto-triplets, but does not divide by the simmetry factor.
"""
function histogram!(xy1, xy2, xy3, dr, hist)

    N = size(hist)[1]
    # Check for identical cells in which case do auto-counts
    for i1 in 1:size(xy1)[2], i2 in 1:size(xy2)[2], i3 in 1:size(xy3)[2]
        p1 = view(xy1, :, i1)
        p2 = view(xy2, :, i2)
        p3 = view(xy3, :, i3)
        r12 = distance(p1, p2)
        r23 = distance(p2, p3)
        r31 = distance(p3, p1)
        if r12 == 0 || r23 == 0 || r31 == 0 continue end
        h1, h2, h3 = hist_index(r12, r23, r31, dr)
        if h1 > N || h2 > N || h3 > N
            continue
        else
            hist[h1, h2, h3] += 1
        end
    end
    
end


"""
    triple_loop!(xy_cube1, xy_cube2, xy_cube3, Ngal1, Ngal2, Ngal3, dr, hist)

    Loop over all triplets of neighbours in the xy_cube grid and histogram them.

    Input:
    xy_cube? - 4D float array of x y coordinates. The first two indeces are for the x and y, the last two are indexing the cell.
    Ngal - 2D Int array of how many galaxies ended up in each cell.
    dr - Float. Bin width.
    hist - 3D Int array of histogrammed separations.

    Looks at all triplets that are separatd by no more than one cell and histograms them.
"""
function triple_loop!(xy_cube1, xy_cube2, xy_cube3, Ngal1, Ngal2, Ngal3, dr, hist)

    Ncube = size(xy_cube1)[end] - 2
    ss = neighbouring_triplets()
    for ix in 2:Ncube+1, iy in 2:Ncube+1
        #println(ix, " ", iy)
        for (jx, jy, kx, ky) in ss
        xy1 = view(xy_cube1, :, 1:Ngal1[ix, iy], ix, iy)
        xy2 = view(xy_cube2, :, 1:Ngal2[ix+jx, iy+jy], ix+jx, iy+jy)
        xy3 = view(xy_cube3, :, 1:Ngal3[ix+kx, iy+ky], ix+kx, iy+ky)
        histogram!(xy1, xy2, xy3, dr, hist)
        end
    end
end


"""
    make_grid(x, y, N)

    Make a grid of x and y coordinates with maximum size of rmax.

    Input:
    x - Float Array. x coordinates.
    y - Float Array. y coordinates.
    N - Int. Grid size.

    Output:
    xy_cube - 4D array of x and y arranged on a grid.
    Ngal - 2D array of integers. Number of galaxies in each gridcell.
    
    The first two indeces in xy_cube are x and y coordinates, the second one goes over galaxies, the last two reference the grid.
    There is a periodic padding around the cubes. The actual cube is between (2,N+1) row/columns 1 is equivalent
    to row/column N+1. Similarly for row/columns N+2 and 2.
"""
function make_cube(x, y, N)
    Ngal = zeros(Int, N+2, N+2)
    xmin = minimum(x)
    xmax = maximum(x)
    ymin = minimum(y)
    ymax = maximum(y)
    L = xmax - xmin

    for (xx, yy) in zip(x, y)
        # 1e-6 to avoid getting zero when xx or yy is exactly at the lower edge.
        Nx = ceil(Int, (xx - xmin + 1e-6)/(xmax - xmin + 1e-6)*N)
        Ny = ceil(Int, (yy - ymin + 1e-6)/(ymax - ymin + 1e-6)*N)
        Ngal[Nx+1, Ny+1] += 1
    end

    # Create periodic padding around the inner cube
    for i in 2:N+1
        Ngal[1,i] = Ngal[N+1,i]
        Ngal[N+2,i] = Ngal[2,i]
        Ngal[i,1] = Ngal[i,N+1]
        Ngal[i,N+2] = Ngal[i,2]
    end
    Ngal[1,1] = Ngal[N+1,N+1]
    Ngal[1,N+2] = Ngal[N+1,2]
    Ngal[N+2,1] = Ngal[2,N+1]
    Ngal[N+2,N+2] = Ngal[2,2]
    # end periodic padding

    Ncounter = copy(Ngal)
    # Cell with the largest number of galaxies
    Nmax = maximum(Ngal)
    # x-index, y-index, gal-index, xy-index
    # This cube is not going to be filled all the way in all cells
    xy_cube = zeros(2, Nmax, N+2, N+2)
    for (xx, yy) in zip(x, y)
        Nx = ceil(Int, (xx - xmin + 1e-6)/(xmax - xmin + 1e-6)*N)
        Ny = ceil(Int, (yy - ymin + 1e-6)/(ymax - ymin + 1e-6)*N)
        xy_cube[:, Ncounter[Nx+1,Ny+1], Nx+1, Ny+1] = [xx, yy]
        Ncounter[Nx+1, Ny+1] -= 1
    end

    # Create periodic padding around the inner cube
    for i in 2:N+1
        xy_cube[:,:,1,i] = xy_cube[:,:,N+1,i] .+ [-L; 0]
        xy_cube[:,:,N+2,i] = xy_cube[:,:,2,i] .+ [L; 0]
        xy_cube[:,:,i,1] = xy_cube[:,:,i,N+1] .+ [0; -L]
        xy_cube[:,:,i,N+2] = xy_cube[:,:,i,2] .+ [0; L]
    end
    xy_cube[:,:,1,1] = xy_cube[:,:,N+1,N+1] .+ [-L; -L]
    xy_cube[:,:,1,N+2] = xy_cube[:,:,N+1,2] .+ [-L; L]
    xy_cube[:,:,N+2,1] = xy_cube[:,:,2,N+1] .+ [L; -L]
    xy_cube[:,:,N+2,N+2] = xy_cube[:,:,2,2] .+ [L; L]
    # end periodic padding

    return xy_cube, Ngal
end


"""
    reduce_hist(hist, dr)

    Input:
    hist - NxNxN array of floats
    dr - bin size (rmin = 0 assumed)

    Output:
    hist_reduced - histogram reduced to a more convenient format ready to print.

    The columns are r1, r2, r3, DDD.
"""
function reduce_hist(hist, dr)
    N = size(hist)[1]
    rmid = collect(range(dr/2, length=N, step=dr))

    Ncol = 0
    for i1 in 1:N, i2 in i1:N, i3 in i2:N
        r1 = rmid[i1] # The shortest
        r2 = rmid[i2]
        r3 = rmid[i3] # The longest
        if r3 + dr/2 < r2 - dr/2 + r1 - dr/2
            Ncol += 1
        end
    end

    hist_reduced = zeros(Ncol, 4)
    Ncol = 1
    for i1 in 1:N, i2 in i1:N, i3 in i2:N
        r1 = rmid[i1] # The shortest
        r2 = rmid[i2]
        r3 = rmid[i3] # The longest
        if r3 + dr/2 < r2 - dr/2 + r1 - dr/2
            hist_reduced[Ncol, 1:3] = [r1 r2 r3]
            if i1 == i2 == i3 # scalene
                hist_reduced[Ncol, 4] = hist[i1, i2, i3]
            elseif i2 == i3 # isoceles
                hist_reduced[Ncol, 4] = hist[i1,i2,i3] + hist[i2,i1,i3] + hist[i2,i3,i1]
            elseif i1 == i2 # isoceles
                hist_reduced[Ncol, 4] = hist[i1,i2,i3] + hist[i1,i3,i2] + hist[i3,i1,i2]
            else # regular
                hist_reduced[Ncol, 4] = hist[i1,i2,i3] + hist[i1,i3,i2] + hist[i2,i1,i3] + hist[i2,i3,i1] + hist[i3,i1,i2] + hist[i3,i2,i1]
            end
            Ncol += 1
        end
    end
    
    return hist_reduced

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
function DDD(x1, y1, x2, y2, x3, y3, dr, Nbin)
    L = maximum([x1 x2 x3 y1 y2 y3]) - minimum([x1 x2 x3 y1 y2 y3])
    Ncell = floor(Int, L/(dr*Nbin))
    xy1_cube, N1_cube = make_cube(x1, y1, Ncell)
    xy2_cube, N2_cube = make_cube(x2, y2, Ncell)
    xy3_cube, N3_cube = make_cube(x3, y3, Ncell)
    hist = zeros(Nbin, Nbin, Nbin)
    triple_loop!(xy1_cube, xy2_cube, xy3_cube, N1_cube, N2_cube, N3_cube, dr, hist)
    rhist = reduce_hist(hist, dr)
end
