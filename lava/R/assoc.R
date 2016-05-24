normal.threshold <- function(object,p=coef(object),...) {
    M <- moments(object,p=p)
    ord <- ordinal(Model(object))
    K <- attributes(ord)$K
    cK <- c(0,cumsum(K-1))
    breaks.orig <- list()
    for (i in seq(K)) {
        breaks.orig <- c(breaks.orig,list(M$e[seq(K[i]-1)+cK[i]]))
    }
    breaks <- lapply(breaks.orig, ordreg_threshold)
    names(breaks) <- names(K)
    ii <- match(names(K),vars(object))
    sigma <- M$Cfull[ii,ii]
    list(breaks=breaks,sigma=sigma,mean=M$v[ii],K=K)
}

prob.normal <- function(sigma,breaks,breaks2=breaks) {
    if (ncol(sigma)!=2 || missing(breaks)) stop("Wrong input")
    P <- matrix(ncol=length(breaks2)-1, nrow=length(breaks)-1)
    for (i in seq(length(breaks)-1))
        for (j in seq(length(breaks2)-1))
            P[i,j] <- mets::pmvn(lower=c(breaks[i],breaks2[j]),upper=c(breaks[i+1],breaks2[j+1]),sigma=sigma)
    return(P)
}

assoc <- function(P,sigma,breaks,...) {
    if (missing(P)) P <- prob.normal(sigma,breaks,...)
    Agree <- sum(diag(P))
    marg.row <- rowSums(P)
    marg.col <- colSums(P)
    Chance <- sum(marg.row*marg.col)
    kap <- (Agree-Chance)/(1-Chance)
    gam <- goodmankruskal_gamma(P)$gamma
    inf <- information_assoc(P)
    res <- c(list(kappa=kap,gamma=gam),inf)
    if (!missing(sigma)) res <- c(res,rho=sigma[1,2])
    return(res)    
}
