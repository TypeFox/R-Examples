imwr.imwd <-
function (imwd, bc = imwd$bc, verbose = FALSE, ...) 
{
    if (verbose == TRUE) 
        cat("Argument checking...")
    ctmp <- class(imwd)
    if (is.null(ctmp)) 
        stop("imwd has no class")
    else if (ctmp != "imwd") 
        stop("imwd is not of class imwd")
    if (imwd$type == "station") 
        stop("Cannot invert nonodecimated wavelet transform using imwr")
    filter <- imwd$filter
    if (verbose == TRUE) 
        cat("...done\nFirst/last database...")
    fl.dbase <- imwd$fl.dbase
    first.last.c <- fl.dbase$first.last.c
    first.last.d <- fl.dbase$first.last.d
    if (verbose == TRUE) 
        cat("...extracted\n")
    ImCC <- imwd$w0Lconstant
    if (verbose == TRUE) 
        cat("Reconstructing...")
    for (level in seq(2, 1 + imwd$nlevels)) {
        if (verbose == TRUE) 
            cat(level - 1, " ")
        LengthCin <- first.last.c[level - 1, 2] - first.last.c[level - 
            1, 1] + 1
        LengthCout <- first.last.c[level, 2] - first.last.c[level, 
            1] + 1
        LengthDin <- first.last.d[level - 1, 2] - first.last.d[level - 
            1, 1] + 1
        error <- 0
        ImOut <- rep(0, LengthCout^2)
        nbc <- switch(bc, periodic = 1, symmetric = 2)
        if (is.null(nbc)) 
            stop("Unknown boundary handling")
        z <- .C("StoIRS", ImCC = as.double(ImCC), ImCD = as.double(imwd[[lt.to.name(level - 
            2, "CD")]]), ImDC = as.double(imwd[[lt.to.name(level - 
            2, "DC")]]), ImDD = as.double(imwd[[lt.to.name(level - 
            2, "DD")]]), LengthCin = as.integer(LengthCin), firstCin = as.integer(first.last.c[level - 
            1, 1]), LengthDin = as.integer(LengthDin), firstDin = as.integer(first.last.d[level - 
            1, 1]), H = as.double(filter$H), LengthH = as.integer(length(filter$H)), 
            LengthCout = as.integer(LengthCout), firstCout = as.integer(first.last.c[level, 
                1]), lastCout = as.integer(first.last.c[level, 
                2]), ImOut = as.double(ImOut), nbc = as.integer(nbc), 
            error = as.integer(error), PACKAGE = "LS2W")
        error <- z$error
        if (error != 0) {
            cat("Error was ", error, "\n")
            stop("Error reported")
        }
        ImCC <- z$ImOut
    }
    if (verbose == TRUE) 
        cat("\nReturning image\n")
    matrix(ImCC, nrow = 2^(imwd$nlevels))
}
