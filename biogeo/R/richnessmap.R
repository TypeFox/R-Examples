richnessmap <-
function (dat, rst, option = "richness") 
{
    f1 <- ifelse(option == "richness", 1, 0)
    f2 <- ifelse(option == "records", 1, 0)
    if ((f1 + f2) < 1) {
        stop("invalid option")
    }
    rst0 <- rst * 0
    v <- values(rst0)
    xy <- data.frame(dat$x, dat$y)
    vals <- extract(rst, xy)
    f <- which(!is.na(vals))
    ce1 <- cellFromXY(rst, xy[f, ])
    dat2 <- data.frame(Species = dat$Species[f], cell = ce1, 
        n = 1)
    if (option == "richness") {
        u <- unique(dat2$cell)
        d <- rep(0, length(u))
        for (i in 1:length(u)) {
            bi <- dat2[dat2$cell == u[i], ]
            d[i] <- length(unique(bi$Species))
        }
        v[u] <- d
        rich <- setValues(rst, v)
    }
    if (option == "records") {
        nrec <- aggregate(n ~ cell, sum, data = dat2)
        v[nrec$cell] <- nrec$n
        rich <- setValues(rst, v)
    }
    return(rich)
}
