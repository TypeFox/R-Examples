findpeaks <- function(x,nups = 1, ndowns = nups, zero = "0", peakpat = NULL, 
                      # peakpat = "[+]{2,}[0]*[-]{2,}", 
                      minpeakheight = -Inf, minpeakdistance = 1,
                      threshold = 0, npeaks = 0, sortstr = FALSE)
{
    stopifnot(is.vector(x, mode="numeric"))
    if (! zero %in% c('0', '+', '-'))
        stop("Argument 'zero' can only be '0', '+', or '-'.")

    # transform x into a "+-+...-+-" character string
    xc <- paste(as.character(sign(diff(x))), collapse="")
    xc <- gsub("1", "+", gsub("-1", "-", xc))
    # transform '0' to zero
    if (zero != '0') xc <- gsub("0", zero, xc)

    # generate the peak pattern with no of ups and downs
    if (is.null(peakpat)) {
        peakpat <- sprintf("[+]{%d,}[-]{%d,}", nups, ndowns)
    }

    # generate and apply the peak pattern
    rc <- gregexpr(peakpat, xc)[[1]]
    if (rc[1] < 0) return(NULL)

    # get indices from regular expression parser
    x1 <- rc
    x2 <- rc + attr(rc, "match.length")
    attributes(x1) <- NULL
    attributes(x2) <- NULL

    # find index positions and maximum values
    n <- length(x1)
    xv <- xp <- numeric(n)
    for (i in 1:n) {
        xp[i] <- which.max(x[x1[i]:x2[i]]) + x1[i] - 1
        xv[i] <- x[xp[i]]
    }

    # eliminate peaks that are too low
    inds <- which(xv >= minpeakheight & xv - pmax(x[x1], x[x2]) >= threshold)

    # combine into a matrix format
    X <- cbind(xv[inds], xp[inds], x1[inds], x2[inds])

    # eliminate peaks that are near by
    if (minpeakdistance < 1)
        warning("Handling 'minpeakdistance < 1' is logically not possible.")

    # sort according to peak height
    if (sortstr || minpeakdistance > 1) {
        sl <- sort.list(X[, 1], na.last = NA, decreasing = TRUE)
        X <- X[sl, , drop = FALSE]
    }

    # return NULL if no peaks
    if (length(X) == 0) return(c())

	# find peaks sufficiently distant
    if (minpeakdistance > 1) {
        no_peaks <- nrow(X)
        badpeaks <- rep(FALSE, no_peaks)

        # eliminate peaks that are close to bigger peaks
        for (i in 1:no_peaks) {
            ipos <- X[i, 2]
            if (!badpeaks[i]) {
                dpos <- abs(ipos - X[, 2])
                badpeaks <- badpeaks | (dpos > 0 & dpos < minpeakdistance)
            }
        }
        # select the good peaks
        X <- X[!badpeaks, ]
    }

    # Return only the first 'npeaks' peaks
    if (npeaks > 0 && npeaks < nrow(X)) {
        X <- X[1:npeaks, , drop = FALSE]
    }

    return(X)
}
