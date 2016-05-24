randomProfile = function(Freqs){
    nLoci = length(Freqs$loci)
    profile = rep(0, 2*nLoci)

    for(nLoc in 1:nLoci){
        i1 = 2*nLoc - 1
        i2 = i1 + 1
        f = Freqs$freqs[[nLoc]]
        profile[i1:i2] = sample(1:length(f), 2, replace = T, prob = f)

        if(profile[i1] > profile[i2]){
            swap = profile[i1]
            profile[i1] = profile[i2]
            profile[i2] = swap
        }
    }

    class(profile) = "profile"
    return(profile)

}
