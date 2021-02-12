using Test

include("../src/wp3.jl")

@test distance(1.0, 1.0, 4.0, 5.0, 2000.0) == 5.0
@test distance(1.0, 1.0, 2000.0 - 2.0, 5.0, 2000.0) == 5.0
@test distance(1.0, 1.0, 4.0, 2000.0 - 3.0, 2000.0) == 5.0
@test distance(1.0, 1.0, 2000.0 - 2.0, 2000.0 - 3.0, 2000.0) == 5.0
@test distance(4.0, 5.0, 1.0, 1.0, 2000.0) == 5.0
@test distance(2000.0 - 2.0, 5.0, 1.0, 1.0, 2000.0) == 5.0
@test distance(4.0, 2000.0 - 3.0, 1.0, 1.0, 2000.0) == 5.0
@test distance(2000.0 - 2.0, 2000.0 - 3.0, 1.0, 1.0, 2000.0) == 5.0

@test hist_index(4, 5, 6, 0.1) == (40, 50, 60)
@test hist_index(0, 4.2, 3.1, 0.2) == (0, 21, 16)

hist = zeros(5, 5, 5)
x = []
y = []
w = []
ind1 = []
ind2 = []
ind3 = []
histogram!(, , , , , , 2000.0, 1.0, hist)
