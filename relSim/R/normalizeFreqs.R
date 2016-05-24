normalizeFreqs = function(Freqs){
    nLoci = length(Freqs$loci)

    for(i in 1:nLoci){
        Freqs$freqs[[i]] = Freqs$freqs[[i]] / sum(Freqs$freqs[[i]])
    }

    return(Freqs)
}
