function RR(rbinedges)
    Nbins = length(rbinedges) - 1
    hist_RR = zeros(Nbins, 2)
    for i in 1:Nbins
        rmid = (rbinedges[i] + rbinedges[i+1])/2.0
        RRvol = 2*π*(rbinedges[i+1]^2 - rbinedges[i]^2)
        hist_RR[i, :] = [rmid, RRvol]
    end
    return hist_RR
end