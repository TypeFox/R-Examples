randomChild = function(profile, Freqs){
    nLoci = length(Freqs$loci)
    profChild = rep(0, 2*nLoci)

    for(nLoc in 1:nLoci){
        f = Freqs$freqs[[nLoc]]
        a = sample(1:length(f), 1, prob = f)
        u = runif(1)
        i1 = 2*nLoc - 1
        i2 = i1 + 1

        if(u < 0.5){
            profChild[i1:i2] = c(profile[i1], a)
        }else{
            profChild[i1:i2] = c(a, profile[i2])
        }

        if(profChild[i1] > profChild[i2]){
            swap = profChild[i1]
            profChild[i1] = profChild[i2]
            profChild[i2] = swap
        }
    }

    class(profChild) = "profile"

    return(profChild)
}
