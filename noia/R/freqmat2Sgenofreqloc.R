freqmat2Sgenofreqloc <-
function (reference = "F2", i = NULL, freqmat = NULL, sinv = TRUE, 
    tol = 1e-08) 
{
    ans <- NULL
    if (reference == "F2") {
        ans$smat <- Sloc("F2")
        ans$genofreq <- matrix(c(0.25, 0.5, 0.25), ncol = 3)
    }
    else if (reference == "UWR") {
        ans$smat <- Sloc("UWR")
        ans$genofreq <- matrix(c(1/3, 1/3, 1/3), ncol = 3)
    }
    else if (reference == "P1") {
        ans$smat <- Sloc("P1")
        ans$genofreq <- matrix(c(1, 0, 0), ncol = 3)
    }
    else if (reference == "P2") {
        ans$smat <- Sloc("P2")
        ans$genofreq <- matrix(c(0, 0, 1), ncol = 3)
    }
    else if (reference == "F1") {
        ans$smat <- Sloc("F1")
        ans$genofreq <- matrix(c(0, 1, 0), ncol = 3)
    }
    else if (reference == "Finf") {
        ans$smat <- Sloc("Finf")
        ans$genofreq <- matrix(c(0.5, 0, 0.5), ncol = 3)
    }
    else if (reference == "noia") {
        if (is.null(freqmat) || is.null(i)) {
            stop("freqmat and i needed for noia model")
        }
        f <- freqmat[i, ]
        if (abs(sum(f) - 1) > tol) {
            stop(paste("frequencies for locus", as.character(i), 
                "do not add up to 1 ([", as.character(f[1]), 
                as.character(f[2]), as.character(f[3]), "])"))
        }
        ans$smat <- matrix(c(1, -f[2] - 2 * f[3], (-2 * f[2] * 
            f[3])/(f[1] + f[3] - (f[1] - f[3])^2), 1, 1 - f[2] - 
            2 * f[3], (4 * f[1] * f[3])/(f[1] + f[3] - (f[1] - 
            f[3])^2), 1, 2 - f[2] - 2 * f[3], (-2 * f[1] * f[2])/(f[1] + 
            f[3] - (f[1] - f[3])^2)), byrow = TRUE, ncol = 3)
        ans$genofreq <- matrix(freqmat[i, ], ncol = 3)
    }
    else if (reference == "G2A") {
        if (is.null(freqmat) || is.null(i)) {
            stop("freqmat and i needed for G2A model")
        }
        p <- freqmat[i]
        ans$smat <- matrix(c(1, -2 * (1 - p), -2 * (1 - p) * 
            (1 - p), 1, -1 + 2 * p, 2 * p * (1 - p), 1, 2 * p, 
            -2 * p * p), byrow = TRUE, ncol = 3)
        ans$genofreq <- matrix(c(p^2, 2 * p * (1 - p), (1 - p)^2), 
            ncol = 3)
    }
    else {
        stop("unknown reference ", reference, ".")
    }
    if (sinv) {
        ans$sinv <- solve(ans$smat)
    }
    return(ans)
}
