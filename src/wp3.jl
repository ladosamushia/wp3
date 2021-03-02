"""
    distance(p1, p2, L)

    Cartesian distance. L-periodic in both dimentions.

    Input:
    p1, p2 - Float{2} contains x and y components.
    L - Float.

    Output:
    dist - Float.
"""
function distance(p1, p2, L)

    dx = abs(p1[1] - p2[1])
    dx > L/2 ? dx -= L : dx = dx
    dy = abs(p1[2] - p2[2])
    dy > L/2 ? dy -= L : dy = dy
    dist = sqrt(dx.^2 + dy.^2)

end

"""
    hist_index(r1, r2, r3, dr)

    Which bin does this r1, r2, r3 combination belongs to if the bin width is dr?
    Inputs:
    all Floats.

    Output:
    triplet of integers.

    Assumes binning is uniform and starts at 0.
"""
function hist_index(r1, r2, r3, dr)

    ind1 = ceil(Int, r1/dr)
    ind2 = ceil(Int, r2/dr)
    ind3 = ceil(Int, r3/dr)
    return ind1, ind2, ind3

end

"""
    histogram!(xy1, xy2, xy3, L, dr, hist)

    look at all triplets in 2D arrays xy1, xy2, xy3.

    Input:
    xy? - Float 2D array.
    L - Float. Box size.
    dr - Float. Bin width.
    hist - Int 3D array.

    xy? can be identical, in which case this computes auto-triplets.
    If only two cells coincide they have to be xy1 and xy2.
"""
function histogram!(xy1, xy2, xy3, L, dr, hist)

    N = size(hist)[1]
    # Check for identical cells in which case do auto-counts
    for i1 in 1:size(xy1)[1]
        p1 = xy1[i1,:]
        if xy1 == xy2
            i2min = i1 + 1
        else
            i2min = 1
        end
        for i2 in i2min:size(xy2)[1]
            p2 = xy2[i2,:]
            if xy2 == xy3
                i3min = i2 + 1
            else
                i3min = 1
            end
            for i3 in i3min:size(xy3)[1]
                p3 = xy3[i3,:]
                r12 = distance(p1, p2, L)
                r23 = distance(p2, p3, L)
                r31 = distance(p3, p1, L)
                h1, h2, h3 = hist_index(r12, r23, r31, dr)
                if h1 > N || h2 > N || h3 > N || h1 < 1 || h2 < 1 || h3 < 1 
                    continue
                else
                    hist[h1, h2, h3] += 1
                end
            end
        end
    end
    
end


"""
    triple_loop!    (xy_cube, Ngal, dr, hist)

    Loop over all triplets of neighbours in the xy_cube grid and histogram them.

    Input:
    xy_cube - 4D float array of x y coordinates. The first two indeces are for the cell, the last two are for x and y
    Ngal - 2D Int array of how many galaxies ended up in each cell.
    dr - Float. Bin width.
    hist - 3D Int array of histogrammed separations.

    Looks at all triplets that are separatd by no more than one cell and histograms them.
"""

function triple_loop!(xy_cube, Ngal, dr, hist)

    step = [0 0 0 0; 0 0 0 1; 0 0 1 -1; 0 0 1 0; 0 0 1 1; 0 1 0 1; 0 1 1 0; 0 1 1 1; 1 -1 1 -1; 1 -1 1 0; 1 0 1 0; 1 0 1 1; 1 1 1 1]

    L = maximum(xy_cube) - minimum(xy_cube)
    Ncube = size(xy_cube)[1]
    for ix in 1:Ncube, iy in 1:Ncube, ss in eachrow(step)
        jx = ix + ss[1] 
        if jx == Ncube + 1 jx = 1 end 
        if jx == 0 jx = Ncube end
        jy = iy + ss[2]        
        if jy == Ncube + 1 jy = 1 end 
        if jy == 0 jy = Ncube end
        kx = ix + ss[3]
        if kx == Ncube + 1 kx = 1 end 
        if kx == 0 kx = Ncube end
        ky = iy + ss[4]
        if ky == Ncube + 1 ky = 1 end 
        if ky == 0 ky = Ncube end

        xy1 = view(xy_cube, ix, iy, 1:Ngal[ix, iy], :)
        xy2 = view(xy_cube, jx, jy, 1:Ngal[jx, jy], :)
        xy3 = view(xy_cube, kx, ky, 1:Ngal[kx, ky], :)
        histogram!(xy1, xy2, xy3, L, dr, hist)
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
    
    The first two indeces in xy_cube reference the grid, the third goes over galaxies, the fourth is x/y.
"""
function make_cube(x, y, N)
    Ngal = zeros(Int, N, N)
    xmin = minimum(x)
    xmax = maximum(x)
    ymin = minimum(y)
    ymax = maximum(y)
    for (xx, yy) in zip(x, y)
        # 1e-6 to avoid getting zero when xx or yy is exactly at the lower edge.
        Nx = ceil(Int, (xx - xmin + 1e-6)/(xmax - xmin + 1e-6)*N)
        Ny = ceil(Int, (yy - ymin + 1e-6)/(ymax - ymin + 1e-6)*N)
        Ngal[Nx, Ny] += 1
    end
    Ncounter = copy(Ngal)
    # Cell with the largest number of galaxies
    Nmax = maximum(Ngal)
    # x-index, y-index, gal-index, xy-index
    # This cube is not going to be filled all the way in all cells
    xy_cube = zeros(N, N, Nmax, 2)
    for (xx, yy) in zip(x, y)
        Nx = ceil(Int, (xx - xmin + 1e-6)/(xmax - xmin + 1e-6)*N)
        Ny = ceil(Int, (yy - ymin + 1e-6)/(ymax - ymin + 1e-6)*N)
        xy_cube[Nx, Ny, Ncounter[Nx, Ny], :] = [xx, yy]
        Ncounter[Nx, Ny] -= 1
    end
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
        if r3 < r2 + r1
            Ncol += 1
        end
    end

    hist_reduced = zeros(Ncol, 4)
    Ncol = 1
    for i1 in 1:N, i2 in i1:N, i3 in i2:N
        r1 = rmid[i1] # The shortest
        r2 = rmid[i2]
        r3 = rmid[i3] # The longest
        if r3 < r2 + r1
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
function DDD(x, y, dr, Nbin)
    L = maximum([x y]) - minimum([x y])
    Ncell = floor(Int, L/(dr*Nbin))
    xy_cube, N_cube = make_cube(x, y, Ncell)
    hist = zeros(Nbin, Nbin, Nbin)
    triple_loop!(xy_cube, N_cube, dr, hist)
    rhist = reduce_hist(hist, dr)
end

function V2cap(r1, r2, r3)
    π*(r2+r1-r3)^2*(r3^2+2*r3*r1-3*r1^2+2*r3*r2+6*r1*r2-3*r2^2)/12/d
end

function RRR(r1, r2, r3, dr)
    r1min = r1 - dr/2
    r1max = r1 + dr/2
    r2min = r2 - dr/2
    r2max = r2 + dr/2
    Vtot = V2cap(r1max, r2max, r3)
    Vtot -= V2cap(r1min, r2max, r3)
    Vtot -= V2cap(r1max, r2min, r3)
    Vtot += V2cap(r1min, r2min, r3)
end