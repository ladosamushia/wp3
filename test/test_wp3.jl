using Test

include("../src/wp3.jl")

@test distance([1.0, 1.0], [4.0, 5.0], 2000.0) == 5.0
@test distance([1.0, 1.0], [2000.0 - 2.0, 5.0], 2000.0) == 5.0
@test distance([1.0, 1.0], [4.0, 2000.0 - 3.0], 2000.0) == 5.0
@test distance([1.0, 1.0], [2000.0 - 2.0, 2000.0 - 3.0], 2000.0) == 5.0
@test distance([4.0, 5.0], [1.0, 1.0], 2000.0) == 5.0
@test distance([2000.0 - 2.0, 5.0], [1.0, 1.0], 2000.0) == 5.0
@test distance([4.0, 2000.0 - 3.0], [1.0, 1.0], 2000.0) == 5.0
@test distance([2000.0 - 2.0, 2000.0 - 3.0], [1.0, 1.0], 2000.0) == 5.0

@test hist_index(4, 5, 6, 0.1) == (40, 50, 60)
@test hist_index(0, 4.2, 3.1, 0.2) == (0, 21, 16)

"""
hist = zeros(5, 5, 5)
x = [1, 2, 3, 4, 5, 6, 7, 8, 1999]
y = [1, 2, 3, 4, 5, 6, 7, 8, 1999]
w = [1, 1, 3, 2, 1, 1, 1, 1, 1]
ind1 = [1, 3, 6]
ind2 = [2, 4, 7]
ind3 = [5, 8, 9]
histogram!(ind1, ind2, ind3, x, y, w, 2000.0, 1.0, hist)
"""

x = [0, 5.1, 10]
y = [0, 4.2, 10]
xyc, Ngal = make_cube(x, y, 10)
for i in 1:10, j in 1:10
    if i == 6 && j == 5
        @test(xyc[i,j,1,:] == [5.1,4.2])
    elseif i == 1 && j == 1
        @test(xyc[i,j,1,:] == [0,0])
    elseif i == 10 && j == 10
        @test(xyc[i,j,1,:] == [10,10])
    else
        @test(xyc[i,j,1,:] == [0,0])
    end
end
@test sum(Ngal) == 3
