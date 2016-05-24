plot.tos <-
function (x, mctype = "FDR", sub = NULL,
	xlab = "Time", arrow.length = 0.05, 
    verbose = FALSE, ...) 
{
    object <- x
    x <- x$x
	if (is.null(sub))
		sub <- paste("MC Type: ", mctype, ". ", object$nreject, " rejected.")
    ts.plot(x, xlab=xlab, sub=sub, col="gray80",...)
    st <- summary.tos(object, mctype = mctype, quiet = TRUE)
    nreject <- st$nreject
    if (nreject == 0) 
        return(NULL)
    st <- st$rejlist
    stPlevs <- stHlevs <- NULL
    for (i in 1:length(st)) {
        stPlevs <- c(stPlevs, st[[i]][1])
        stHlevs <- c(stHlevs, st[[i]][2])
    }
    nPlevs <- length(min(stPlevs):max(stPlevs))
    nHlevs <- length(min(stHlevs):max(stHlevs))
    ry <- range(x)
    mny <- ry[1]
    mxy <- ry[2]
    mainy <- seq(from = mny, to = mxy, length = nPlevs + 1)
    lyn <- min(stHlevs)
    lyx <- max(stHlevs)
    littley <- seq(from = 0, to = (mainy[2] - mainy[1]), length = lyx - 
        lyn + 2)
    if (verbose == TRUE) {
        cat("nPlevs: ", nPlevs, "\n")
        cat("mny, mxy: ", mny, mxy, "\n")
        cat("mainy: ")
        print(mainy)
        cat("littley: ")
        print(littley)
    }
    abline(h = mainy[1:(length(mainy) - 1)], lty = 2)
    axis(4, at = mainy[1:(length(mainy) - 1)], labels = min(stPlevs):max(stPlevs))
    J <- IsPowerOfTwo(length(x))
    for (i in 1:length(st)) {
        stP <- st[[i]][1]
        stH <- st[[i]][2]
        ix <- st[[i]][c(-1, -2)]
        for (j in 1:length(ix)) {
            xl <- 2^(J - stH) * (ix[j] - 1)
            xr <- 2^(J - stH) * (ix[j])
            yy <- mainy[stP - min(stPlevs) + 1] + littley[stH - 
                lyn + 1]
            arrows(x0 = xl, x1 = xr, y0 = yy, y1 = yy, code = 3, 
                col = 2, length = arrow.length)
            if (verbose == TRUE) {
                cat("stP, stH: ", stP, stH, "\n")
                cat("[xl, xt] ", xl, xr, mainy[stP - min(stPlevs) + 
                  1], littley[stH - lyn + 1], "\n")
                scan()
            }
        }
    }
}
