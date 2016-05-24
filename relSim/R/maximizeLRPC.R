maximizeLRPC = function(N, Freqs,
                        P1 = randomProfile(Freqs),
                        P2 = randomProfile(Freqs), BlockSize = N/10){

    f = unlist(Freqs$freqs)
    alleleCounts = sapply(Freqs$freqs, length)
    nLoci = length(Freqs$loci)

    lrMax = lrPC(P1, P2, Freqs)

    nBlocks = N/BlockSize

    pb = txtProgressBar(min = 0, max = nBlocks, style = 3)

    for(block in 1:nBlocks){
        u = runif(4*nLoci*BlockSize)

        res = .C("maximizeLRPC", P1 = as.integer(P1), P2 = as.integer(P2),
                 nLoc = as.integer(nLoci),
                 f = as.double(f),
                 nA = as.integer(alleleCounts),
                 u = as.double(u),
                 bs = as.integer(BlockSize), lrMax = as.double(lrMax))

        P1 = res$P1
        class(P1) = "profile"
        P2 = res$P2
        class(P2) = "profile"
        lrMax = res$lrMax

        #print(P1, horiz = T)
        #print(P2, horiz = T)
        #cat(paste(lrMax,"\n"))

        setTxtProgressBar(pb, block)
    }

    return(list(P1 = P1, P2 = P2, lrMax = lrMax))
}
