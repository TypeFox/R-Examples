speciescount <-
function (dat, orderby = "Species") 
{
    spp <- unique(dat$Species)
    nsp <- {
    }
    for (i in 1:length(spp)) {
        fsp <- which(dat$Species == spp[i])
        sp <- dat[fsp, ]
        n <- nrow(sp)
        excl <- sum(sp$Exclude)
        nu <- n - excl
        d <- data.frame(Species = spp[i], ntot = n, nuniq = nu)
        nsp <- rbind(nsp, d)
    }
    if (orderby == "ntot") {
        f <- order(nsp$ntot, decreasing = TRUE)
        nsp <- nsp[f, ]
    }
    if (orderby == "nuniq") {
        f <- order(nsp$nuniq, decreasing = TRUE)
        nsp <- nsp[f, ]
    }
    if (orderby == "Species") {
        f <- order(spp)
        nsp <- nsp[f, ]
    }
    return(nsp)
}
