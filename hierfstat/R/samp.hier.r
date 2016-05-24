"samp.between" <-
function (lev) 
{
    y <- 1:length(lev)
    nlev <- nlevels(factor(lev))
    nl <- 1:nlev
    x <- list()
    for (i in 1:nlev) x[[i]] <- y[lev == i]
    dum <- sample(nl)
    return(unlist(x[dum]))
}
"samp.between.within" <-
function (inner.lev, outer.lev) 
{
    y <- 1:length(inner.lev)
    nlev.o <- nlevels(factor(outer.lev))
    z <- NULL
    for (j in 1:nlev.o) {
        lev.i <- as.integer(levels(factor(inner.lev[outer.lev == 
            j])))
        x <- list()
        if (length(lev.i) == 1) {
            z[[j]] <- y[inner.lev == lev.i]
        }
        else {
            for (i in 1:length(lev.i)) x[[i]] <- y[inner.lev == 
                lev.i[i]]
            dum <- sample(1:length(lev.i))
            z[[j]] <- unlist(x[dum])
        }
    }
    return(unlist(z))
}
"samp.within" <-
function (lev) 
{
    y <- 1:length(lev)
    nlev <- nlevels(factor(lev))
    nl <- as.integer(factor(lev))
    x <- list()
    for (i in 1:nlev) if (length(y[nl==i])>1) x[[i]] <- sample(y[nl == i]) else x[[i]]<-y[nl==i]
    return(unlist(x))
}
