include("wp3.jl")
include("random_wp3.jl")

function wp3()
    DDD_counts = DDD(x1, y1, x2, y2, x3, y3, dr, Nbin)
    RRR_counts = RRR(rbinedges)
    wp3 = DDD_counts./RRR_counts - 1
end