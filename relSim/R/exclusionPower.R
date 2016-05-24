exclusionPower = function(Freqs){
    ep = 1 - sapply(Freqs$freqs,function(p){2 * sum(p^2)^2 - sum(p^4)})
    names(ep) = Freqs$loci
    return(ep)
}

ep = function(Freqs){
    exclusionPower(Freqs)
}


