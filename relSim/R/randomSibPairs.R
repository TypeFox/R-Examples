randomSibPairs = function(Freqs, BlockSize = 1){
    ## nLoci = length(Freqs$loci)
    ## profSib = matrix(0, nc = 2, nr = nLoci)

    ## for(nLoc in 1:nLoci){
    ##     f = Freqs$freqs[[nLoc]]
    ##     i = sample(1:4, 1)
    ##     a = sample(1:length(f), 2, replace = TRUE, prob = f)

    ##     switch(i,
    ##            {profSib[nLoc,] = profile[nLoc,]},
    ##            {profSib[nLoc,] = c(profile[nLoc,1], a[1])},
    ##            {profSib[nLoc,] = c(a[1], profile[nLoc,2])},
    ##            {profSib[nLoc,] = a}
    ##            )


    ##     if(profSib[nLoc,1] > profSib[nLoc, 2]){
    ##         swap = profSib[nLoc, 1]
    ##         profSib[nLoc, 1] = profSib[nLoc, 2]
    ##         profSib[nLoc, 2] = swap
    ##     }
    ## }

    nLoci = length(Freqs$loci)
    Profile = vector(mode = "list", length = BlockSize)
    
    Sib1 = randomProfiles(Freqs$freqs, BlockSize);
    Sib2 = randomSibs(Sib1, Freqs$freqs, BlockSize);

    for(b in 1:BlockSize){
        i1 = (b - 1)*2*nLoci + 1
        i2 =  b * 2*nLoci

        Profile[[b]] = vector(mode = "list", length = 2)
        names(Profile[[b]]) = c("sib1", "sib2")

        Profile[[b]]$sib1 = Sib1[i1:i2]
        class(Profile[[b]]$sib1) = "profile"

        Profile[[b]]$sib2 = Sib2[i1:i2]
        class(Profile[[b]]$sib2) = "profile"
    }

    if(BlockSize==1){
        return(Profile[[1]])
    }else{
        return(Profile)
    }
}
