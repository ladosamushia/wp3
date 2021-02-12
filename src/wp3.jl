"""
    distance(x1, x2, y1, y2)

    Cartesian distance.

    Input:
    x1, x2, y1, y2 - Floats.

    Output:
    dist - Float.
"""
function distance(x1, x2, y1, y2)

    return dist
end

"""
    histogram!(ind1, ind2, ind3, x, y, hist)

    look at all triplets in ind1, ind2, ind3 and histogram the distances of x and y that they point to.

    Input:
    ind1 - Int array.
    ind2 - Int array.
    Ind3 - Int array.
    x - Float array.
    y - Float array.
    hist - Float 3D array.

    ind1, ind2, ind3 can be identical, in which case this computes auto-triplets.
"""
function histogram!(x1, x2, x3, y1, y2, y3, hist)

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
