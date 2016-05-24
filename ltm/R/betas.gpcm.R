betas.gpcm <-
function (thetas, nitems, ncatg, constraint, keep.names = FALSE) {
    betas <- if (constraint == "gpcm") {
        ii <- rep(1:nitems, ncatg)
        split(thetas, ii)
    } else if (constraint == "1PL") {
        nt <- length(thetas)
        ii <- rep(1:nitems, ncatg - 1)
        lapply(split(thetas[-nt], ii), function (x) c(x, thetas[nt]))
    } else {
        ii <- rep(1:nitems, ncatg - 1)
        lapply(split(thetas, ii), function (x) c(x, 1))       
    }
    if (!keep.names)
        names(betas) <- NULL
    betas
}
