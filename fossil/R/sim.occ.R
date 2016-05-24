sim.occ <-
function(total.species = 100, endemics = 0.1, regions = 3, locs = 30, avg.abund = 1) {
    cosmop <- round(total.species*(1-endemics))
    endem <- total.species-cosmop
    end.per.reg <- floor(endem/regions)
    extras <- endem%%regions
    orig.groups <- matrix(0, total.species, regions*locs)
    reg<-1
    for (i in 1:(regions*locs)) {
        #cosmopolitan spp
        orig.groups[1:cosmop, i] <- 1
        #endemic spp start at
        strt <- cosmop+1+((reg-1)*end.per.reg)
        ends <- strt+end.per.reg-1
        if (reg==regions) ends <- ends+extras
        orig.groups[strt:ends, i] <- 1
        orig.groups[orig.groups[,i]==1,i]<-round(rlnorm(length(orig.groups[orig.groups[,i]==1,i]), 0, 1)*avg.abund)
        if (i%%locs==0) reg<-reg+1
    }  
    return(orig.groups)
}
