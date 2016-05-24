randomSib = function(profile, Freqs){
    nLoci = length(Freqs$loci)
    profSib = rep(0, 2*nLoci)

    for(nLoc in 1:nLoci){
        f = Freqs$freqs[[nLoc]]
        i = sample(1:4, 1)
        a = sample(1:length(f), 2, replace = TRUE, prob = f)
        i1 = 2 * nLoc - 1
        i2 = i1 + 1

        switch(i,
               {profSib[i1:i2] = profile[i1:i2]},
               {profSib[i1:i2] = c(profile[i1], a[1])},
               {profSib[i1:i2] = c(a[1], profile[i2])},
               {profSib[i1:i2] = a}
               )


        if(profSib[i1] > profSib[i2]){
            swap = profSib[i1]
            profSib[i1] = profSib[i2]
            profSib[i2] = swap
        }
    }
    class(profSib) = "profile"
    return(profSib)
}
