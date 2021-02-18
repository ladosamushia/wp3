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
    for i1 in 1:length(xy1)[2]
        if xy1 == xy2
            i2min = i1 + 1
        else
            i2min = 1
        end
        for i2 in i2min:length(xy2)[2]
            if xy2 == xy3
                i3min = i2 + 1
            else
                i3min = 1
            end
            for i3 in i3min:length(xy3)[2]
                r12 = distance(p1, p2, L)
                r23 = distance(p2, p3, L)
                r31 = distance(p3, p1, L)
                h1, h2, h3 = hist_index(r12, r23, r13, dr)
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
    triple_loop!    (xy_cube, dr, hist)

    Loop over all triplets of neighbours in the xy_cube grid and histogram them.

    Input:
    xy_cube - 2D Int array of indeces.
    dr - Float. Bin width.
    hist - 3D Int array of histogrammed separations.

    Looks at all triplets that are separatd by no more than one cell and histograms them.
"""
function triple_loop!(xy_cube, dr, hist)
    L = maximum(xy_cube) - minimum(xy_cube)
    Nmax = size(hist)[1]
    Ncube = size(xy_cube)[1]
    for i1 in 1:Ncube, j1 in 1:Ncube
        for i2 in 1:Ncube, j2 in 1:Ncube
            di12 = abs(i1 - i2)
            dj12 = abs(j1 - j2)
            # Check these are neighbouring cells
            if (di12 == 1 || di12 == Ncube - 1) && (dj12 == 1 || dj12 == Ncube - 1) 
                for i3 in 1:Ncube, j3 in 1:Ncube
                    di23 = abs(i2 - i3)
                    di13 = abs(i1 - i3)
                    dj23 = abs(j2 - j3)
                    dj13 = abs(j1 - j3)
                    # Check these are neighbouring cells
                    if (di23 == 1 || di23 == Ncube - 1) && (dj23 == 1 || dj23 == Ncube - 1)  && (di13 == 1 || di13 == Ncube - 1) && (dj13 == 1 || dj13 == Ncube - 1)
                        histogram!(xy_cube[i1, j1], xy_cube[i2, j2], xy_cube[i3, j3], L, dr, hist)
                    end
                end
            end
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
