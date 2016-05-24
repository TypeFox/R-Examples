lkappa <-
function(r, type="Cohen", weights="squared") 
{
    nrater <- dim(r)[2]
    kappas <- vector(length = nrater * (nrater - 1)/2)
    c <- 0
    for (i in 2:nrater) for (j in 1:(i - 1))
    {
        c <- c + 1
        if (type == "Cohen") kappas[c] <- ckappa(r[, c(i, j)])[[2]]
        else  kappas[c] <- wkappa(r[, c(i, j)], weights=weights)[[3]]
    }
    return(mean(kappas))
}
