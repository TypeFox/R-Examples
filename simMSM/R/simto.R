simto <- function(entry.ij, from.ij, mpl, eta.ij, x.i, max.time, pme){
    exit.ij <- simexit(entry.ij, all.bhr = mpl[[from.ij]]$bhr, x.i = x.i,
                       eta.ij = eta.ij, max.time = max.time, pme = pme)$new.exit
    hr.at.exit.ij <- rep(NA, length(eta.ij))
    for(hi in mpl[[from.ij]]$all.to){
        hr.at.exit.ij[hi] <- hr(bhr = mpl[[from.ij]]$bhr[[hi]],
                                t = exit.ij, 
                                eta.ij = eta.ij[[hi]], 
                                x.i = x.i) * pme[hi]
    }
    hr.at.exit.ij <- hr.at.exit.ij[!is.na(hr.at.exit.ij)]
    if(length(hr.at.exit.ij) > 1.5){
        probs <- hr.at.exit.ij/sum(hr.at.exit.ij)
        to.ij <- sample(mpl[[from.ij]]$all.to, size = 1, prob = probs)
    }else{
        to.ij <- as.numeric(mpl[[from.ij]]$all.to)
    }
    return(list(entry.ij = entry.ij, 
                exit.ij = exit.ij, 
                from.ij = from.ij,
                to.ij = to.ij))}