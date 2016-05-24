plot.tosANYN <-
function (x, sub = NULL, xlab = "Time", arrow.length = 0.05, 
    verbose = FALSE, ...) 
{
    object <- x
    x <- x$x
    if (is.null(sub)) 
        sub <- paste("MC Type: ", object$mc.method[1], ". ", object$nreject, 
            " rejected.")
    plot(x=1:length(x), y=x, xlab = xlab, sub = sub, col = "gray80", type="l", ...)
    st <- summary.tosANYN(object, quiet=TRUE)
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
#
#   Now have to work out J equivalent here. In plot.tos we used
#   J <- IsPowerOfTwo(length(x))
#
#   This won't work here because length of x is not a power of two. So,
#   replace with a nearby version
#
    lx <- length(x)
    J <- floor(log(lx)/log(2)) 

    if (verbose==TRUE)
	cat("J for non-dyadic length is: ", J, "\n")

#	stP is the scale within the periodogram - the "main" scale
#	stH is the scale within each periodogram scale, the minor scale

    for (i in 1:length(st)) {
        stP <- st[[i]][1]
        stH <- st[[i]][2]
	lv <- st[[i]][3]
        ix <- st[[i]][c(-1, -2, -3)]

	if (verbose==TRUE)
		cat("stP: ", stP, " stH: ", stH, "lv: ", lv, " ix: ", ix, "\n")
        for (j in 1:length(ix)) {

	    mylv <- 2^stH
	
            xl <- lx* (ix[j] - 1)/lv
            xr <- lx*(ix[j])/lv
            yy <- mainy[stP - min(stPlevs) + 1] + littley[stH - 
                lyn + 1]
            arrows(x0 = xl, x1 = xr, y0 = yy, y1 = yy, code = 3, 
                col = 2, length = arrow.length)
            if (verbose == TRUE) {
                cat("stP, stH: ", stP, stH, "\n")
                cat("[xl, xt] ", xl, xr, " mainy: ", mainy[stP - min(stPlevs) + 
                  1], "littley:", littley[stH - lyn + 1],
                " sum: ", mainy[stP - min(stPlevs) + 1] + littley[stH - lyn + 1], "\n")
                scan()
            }
        }
    }
}
