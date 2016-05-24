 Isoqqstat <- function (x, y, fudge, niter) 
{
    xiter.index <- t(sapply(1:niter, function(i) sample(x)))
    to1 <- to2 <- to3 <- to4 <- to5 <- matrix(0, nrow(y), niter)
    fudge.factor <- switch(fudge, pooled = Isofudge(x, y), none = c(rep(0, 
        5)))
    total <- niter
    pb <- tkProgressBar(title = "progress bar", min = 0, max = total, 
        width = 300)
    for (i in 1:niter) {
        yyy0 <- IsoGenemSAM(xiter.index[i, ], as.matrix(y), fudge.factor)
        to1[, i] <- sort(yyy0[[1]])
        to2[, i] <- sort(yyy0[[2]])
        to3[, i] <- sort(yyy0[[3]])
        to4[, i] <- sort(yyy0[[4]])
        to5[, i] <- sort(yyy0[[5]])
        setTkProgressBar(pb, i, title = paste(round(i/total * 
            100, 0), "% done"))
    }
    close(pb)
    L <- IsoGenemSAM(x, as.matrix(y), fudge.factor)
    d <- L[[1]]
    d.sort.list <- sort.list(d)
    d.sort <- d[d.sort.list]
    perm.mean <- rowMeans(to1)
    aa1 = cbind(d.sort, perm.mean, d.sort - perm.mean, d.sort.list)
    row.names(aa1) <- row.names(y)[d.sort.list]
    d <- L[[2]]
    d.sort.list <- sort.list(d)
    d.sort <- d[d.sort.list]
    perm.mean <- rowMeans(to2)
    aa2 = cbind(d.sort, perm.mean, d.sort - perm.mean, d.sort.list)
    row.names(aa2) <- row.names(y)[d.sort.list]
    d <- L[[3]]
    d.sort.list <- sort.list(d)
    d.sort <- d[d.sort.list]
    perm.mean <- rowMeans(to3)
    aa3 = cbind(d.sort, perm.mean, d.sort - perm.mean, d.sort.list)
    row.names(aa3) <- row.names(y)[d.sort.list]
    d <- L[[4]]
    d.sort.list <- sort.list(d)
    d.sort <- d[d.sort.list]
    perm.mean <- rowMeans(to4)
    aa4 = cbind(d.sort, perm.mean, d.sort - perm.mean, d.sort.list)
    row.names(aa4) <- row.names(y)[d.sort.list]
    d <- L[[5]]
    d.sort.list <- sort.list(d)
    d.sort <- d[d.sort.list]
    perm.mean <- rowMeans(to5)
    aa5 <- cbind(d.sort, perm.mean, d.sort - perm.mean, d.sort.list)
    row.names(aa5) <- row.names(y)[d.sort.list]
    res <- list(aa1 = aa1, to1 = to1, aa2 = aa2, to2 = to2, aa3 = aa3, 
        to3 = to3, aa4 = aa4, to4 = to4, aa5 = aa5, to5 = to5)
    return(res)
}
