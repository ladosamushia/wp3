using Test

include("../src/wp3.jl")

# distance

@test distance([1.0, 1.0], [4.0, 5.0], 2000.0) == 5.0
@test distance([1.0, 1.0], [2000.0 - 2.0, 5.0], 2000.0) == 5.0
@test distance([1.0, 1.0], [4.0, 2000.0 - 3.0], 2000.0) == 5.0
@test distance([1.0, 1.0], [2000.0 - 2.0, 2000.0 - 3.0], 2000.0) == 5.0
@test distance([4.0, 5.0], [1.0, 1.0], 2000.0) == 5.0
@test distance([2000.0 - 2.0, 5.0], [1.0, 1.0], 2000.0) == 5.0
@test distance([4.0, 2000.0 - 3.0], [1.0, 1.0], 2000.0) == 5.0
@test distance([2000.0 - 2.0, 2000.0 - 3.0], [1.0, 1.0], 2000.0) == 5.0

# hist_index

@test hist_index(4, 5, 6, 0.1) == (40, 50, 60)
@test hist_index(0, 4.2, 3.1, 0.2) == (0, 21, 16)

# histogram!(xy1, xy2, xy3, L, dr, hist)

hist = zeros(Int, 10, 10, 10)
xy1 = [1 2 3 1999; 1 2 3 1999]
xy2 = [1.01 2.01 3.01 1999.01; 1.01 2.01 3.01 1999.01]
xy3 = [7 8 9; 7 8 9]

histogram!(xy3, xy3, xy3, 2000.0, 1.0, hist)
for i in 1:10, j in 1:10, k in 1:10
    if i == 2 && j == 2 && k == 3
        @test hist[i, j, k] + hist[i, k, j] + hist[j, i, k] + hist[j, k, i] + hist[k, i, j] + hist[k, j, i] == 3
    else
        @test hist[i, j, k] == 0
    end
end

hist = zeros(Int, 10, 10, 10)
histogram!(xy3, xy3, xy1, 2000.0, 1.0, hist)
hist_comp = zeros(Int, 10, 10, 10)
hist_comp[2, 6, 7] = 1
hist_comp[2, 7, 8] = 2
hist_comp[2, 8, 10] = 2
hist_comp[3, 6, 10] = 10
for i in 1:10, j in 1:10, k in 1:10
    hist[i, j, k] = hist_comp[i, j, k]
end

hist = zeros(Int, 10, 10, 10)
histogram!(xy1, xy2, xy3, 2000.0, 1.0, hist)
hist_comp = zeros(Int, 10, 10, 10)
hist_comp[1, 6, 6] = 1
hist_comp[1, 8, 8] = 2
hist_comp[1, 9, 9] = 2
hist_comp[1, 10, 10] = 2
hist_comp[2, 6, 8] = 1
hist_comp[2, 8, 6] = 1
hist_comp[2, 8, 9] = 2
hist_comp[2, 9, 8] = 2
hist_comp[2, 10, 9] = 2
hist_comp[2, 9, 10] = 2
hist_comp[3, 6, 9] = 1
hist_comp[3, 9, 6] = 1
hist_comp[3, 8, 10] = 1
hist_comp[3, 10, 8] = 1
for i in 1:10, j in 1:10, k in 1:10
    hist[i, j, k] = hist_comp[i, j, k]
end

# make_cube

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

# triple_loop
x = [0, 1000, 499, 499, 499, 500, 500, 500, 501, 501, 501]
y = [0, 1000, 499, 500, 501, 499, 500, 501, 499, 500, 501]
hist = zeros(Int, 10, 10, 10)
xy_cube, Ngal = make_cube(x, y, 100)
triple_loop!(xy_cube, Ngal, 1.0, hist)
@test hist[1,1,2] + hist[1,2,1] + hist[2,1,1] == 22
@test sum(hist) == 84

x = [999, 999, 999, 0, 0, 1000, 1, 1, 1]
y = [999, 1000, 1, 999, 0, 1, 999, 1000, 1]
hist = zeros(Int, 10, 10, 10)
xy_cube, Ngal = make_cube(x, y, 100)
triple_loop!(xy_cube, Ngal, 1.0, hist)
@test hist[1,1,2] + hist[1,2,1] + hist[2,1,1] == 22
@test sum(hist) == 84

# reduce_hist
hist = ones(5, 5, 5)
rhist = reduce_hist(hist, 1.0)
for rr in eachrow(rhist)
    if rr[1] == rr[2] == rr[3]
        @test rr[4] == 1
    elseif rr[2] == rr[3] || rr[1] == rr[2]
        @test rr[4] == 3
    else
        @test rr[4] == 6
    end
end

# DDD
x = [0, 1000, 499, 499, 499, 500, 500, 500, 501, 501, 501]
y = [0, 1000, 499, 500, 501, 499, 500, 501, 499, 500, 501]
rhist = DDD(x, y, 1.0, 10)
for rr in eachrow(rhist)
    if rr[1] == 0.5 && rr[2] == 0.5 && rr[3] == 1.5
        @test rr[4] == 22
    end
end
#@test sum(rhist[:,4]) == 84