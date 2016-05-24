outliers <-
function (rid, species, dups, ev) 
{
    uspp <- unique(species)
    nr <- length(species)
    ee <- rep(0, nr)
    ee2 <- rep(0, nr)
    for (j in 1:length(uspp)) {
        spp <- uspp[j]
        fsp <- which(species == spp & dups == 0 & !is.na(ev))
        if (length(fsp) >= 10) {
            xc <- ev[fsp]
            ri <- rid[fsp]
            b1 <- boxplot.stats(xc, coef = 1.5)
            xr <- range(b1$stats)
            fe <- which(xc > xr[2] | xc < xr[1])
            ff1 <- ri[fe]
            fe2 <- rjack(xc)
            ff2 <- ri[fe2]
            ee[ff1] <- 1
            ee2[ff2] <- 1
        }
    }
    out <- cbind(ee, ee2)
    return(out)
}
