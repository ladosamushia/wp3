"""
    distance(x1, y1, x2, y2, L)

    Cartesian distance. L-periodic in both dimentions.

    Input:
    x1, x2, y1, y2 - Floats.
    L - Float.

    Output:
    dist - Float.
"""
function distance(x1, y1, x2, y2, L)

    dx = abs(x2 - x1)
    dx > L/2 ? dx -= L : dx = dx
    dy = abs(y2 - y1)
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
    histogram!(ind1, ind2, ind3, x, y, w, L, dr, hist)

    look at all triplets in ind1, ind2, ind3 and histogram the distances of x and y that they point to.

    Input:
    ind1 - Int array.
    ind2 - Int array.
    Ind3 - Int array.
    x - Float array.
    y - Float array.
    w - Float array. Weights.
    L - Float. Box size.
    dr - Float. Bin width.
    hist - Float 3D array.

    ind1, ind2, ind3 can be identical, in which case this computes auto-triplets.
"""
function histogram!(ind1, ind2, ind3, x, y, w, L, dr, hist)

    N = size(hist)[1]
    for i1 in ind1, i2 in ind2, i3 in ind3
        r12 = distance(x[i1], y[i1], x[i2], y[i2], L)
        r23 = distance(x[i2], y[i2], x[i3], y[i3], L)
        r31 = distance(x[i3], y[i3], x[i1], y[i1], L)
        h1, h2, h3 = hist_index(r12, r23, r13, dr)
        if h1 > N || h2 > N || h3 > N || h1 < 1 || h2 < 1 || h3 < 1 
            continue
        else
            hist[h1, h2, h3] += w[i1]*w[i2]*w[i3]
        end
    end

end

"""
    triple_loop(xy_cube, x, y)

    Loop over all triplets of neighbours in the xy_cube grid and histogram them.

    Input:
    xy_cube - 2D Int array of indeces.
    x - Float array. x coordinates.
    y - Float array. y coordinates.

    Output:
    hist - histogram of triplet separations.

    Looks at all triplets that are separatd by no more than one cell and histograms them.
"""
function triple_loop(xy_cube, x, y)

    return hist
end

"""
    make_grid(x, y, rmax)

    Make a grid of x and y coordinates with maximum size of rmax.

    Input:
    x - Float Array. x coordinates.
    y - Float Array. y coordinates.
    rmax - Float. Maximum separation for the grid.

    Output:
    xy_cube - 2D Int array of indeces.
    
    xy_cube elements are arrays of integers that contain indexes of x and y that fall into that cube.
"""
function make_cube(x, y, rmax)

    return xy_cube
end
