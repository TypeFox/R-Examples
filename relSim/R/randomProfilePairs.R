randomProfilePairs = function(Freqs, BlockSize = 1){
    ## nLoci = length(Freqs$loci)
    ## profile = matrix(0, nc = 2, nr = nLoci)

    ## for(nLoc in 1:nLoci){
    ##     f = Freqs$freqs[[nLoc]]
    ##     profile[nLoc,] = sample(1:length(f), 2, replace = T, prob = f)

    ##     if(profile[nLoc,1] > profile[nLoc, 2]){
    ##         swap = profile[nLoc, 1]
    ##         profile[nLoc, 1] = profile[nLoc, 2]
    ##         profile[nLoc, 2] = swap
    ##     }
    ## }

    ## class(profile) = "profile"
    ## return(profile)

    

    prof1 = randomProfiles(Freqs$freqs, BlockSize)
    prof2 = randomProfiles(Freqs$freqs, BlockSize)
    nLoci = length(Freqs$loci)
    Profile = vector(mode = "list", length = BlockSize)
    
    for(b in 1:BlockSize){
        i1 = (b - 1) * 2 * nLoci + 1
        i2 =  b * 2 * nLoci

        Profile[[b]]$prof1 = prof1[i1:i2]
        Profile[[b]]$prof2 = prof2[i1:i2]
        class(Profile[[b]]$prof1) = "profile"
        class(Profile[[b]]$prof2) = "profile"
    }

    if(BlockSize==1){
        return(Profile[[1]])
    }else{
        return(Profile)
    }
}
