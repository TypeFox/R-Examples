mkcoef <-
function (J, filter.number = 10, family = "DaubLeAsymm") 
{
    if (J >= 0) 
        stop("J must be negative integer")
    coefvec <- NULL
    lvec <- NULL
    BigJ <- whichlevel(-J, filter.number = filter.number, family = family)
    n <- 2^BigJ
    ixvec <- cumsum(2^((BigJ - 1):0))
    tmpwd <- wd(rep(0, n), filter.number = filter.number, family = family)
    for (i in 1:(-J)) {
        tmpwd$D <- rep(0, n - 1)
        tmpwd$D[ixvec[i]] <- 1
        vec <- wr(tmpwd)
        vec <- guyrot(vec, idlastzero(vec))
        vec[abs(vec) < 1e-100] <- 0
        vec <- vec[vec != 0]
        coefvec <- c(coefvec, vec)
        lvec <- c(lvec, length(vec))
    }
    v <- vector("list", -J)
    lh <- 1
    for (i in 1:(-J)) {
        v[[i]] <- coefvec[lh:(lh + lvec[i] - 1)]
        lh <- lh + lvec[i]
    }
    v
}
