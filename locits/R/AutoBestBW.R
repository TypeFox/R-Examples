AutoBestBW <-
function (x, filter.number = 1, family = "DaubExPhase", smooth.dev = var, 
    AutoReflect = TRUE, tol = 0.01, maxits = 200, 
    plot.it = FALSE, verbose = 0, ReturnAll=FALSE) 
{
    EWS <- ewspec3(x = x, filter.number = filter.number, family = family, 
        smooth.dev = smooth.dev,
	AutoReflect = AutoReflect, WPsmooth.type = "wavelet")$S
    J <- EWS$nlevels
    Jmax <- J - 1
    Jmin <- round(J/2)
    specerr <- function(S1, S2, levs) {
        ans <- 0
        for (i in 1:length(levs)) {
            d1 <- accessD(S1, lev = levs[i])
            d2 <- accessD(S2, lev = levs[i])
            ans <- ans + sum((d1 - d2)^2)
        }
    return(ans)
    }
    its <- 0
    R <- 0.61803399
    C <- 1 - R
    ax <- 10
    cx <- round(0.6 * length(x))
    bx <- round(cx/2)
    x0 <- ax
    x3 <- cx
    if (abs(cx - bx) > abs(bx - ax)) {
        x1 <- bx
        x2 <- round(bx + C * (cx - bx))
    }
    else {
        x2 <- bx
        x1 <- round(bx - C * (bx - ax))
    }
    fa.ews <- ewspec3(x = x, filter.number = filter.number, family = family, 
        smooth.dev = smooth.dev, AutoReflect = AutoReflect, WPsmooth.type = "RM", 
        binwidth = ax)$S
    fb.ews <- ewspec3(x = x, filter.number = filter.number, family = family, 
        smooth.dev = smooth.dev, AutoReflect = AutoReflect, WPsmooth.type = "RM", 
        binwidth = bx)$S
    fc.ews <- ewspec3(x = x, filter.number = filter.number, family = family, 
        smooth.dev = smooth.dev, AutoReflect = AutoReflect, WPsmooth.type = "RM", 
        binwidth = cx)$S
    f1.ews <- ewspec3(x = x, filter.number = filter.number, family = family, 
        smooth.dev = smooth.dev, AutoReflect = AutoReflect, WPsmooth.type = "RM", 
        binwidth = x1)$S
    f2.ews <- ewspec3(x = x, filter.number = filter.number, family = family, 
        smooth.dev = smooth.dev, AutoReflect = AutoReflect, WPsmooth.type = "RM", 
        binwidth = x2)$S
    fa <- specerr(fa.ews, EWS, Jmin:Jmax)
    fb <- specerr(fb.ews, EWS, Jmin:Jmax)
    fc <- specerr(fc.ews, EWS, Jmin:Jmax)
    f1 <- specerr(f1.ews, EWS, Jmin:Jmax)
    f2 <- specerr(f2.ews, EWS, Jmin:Jmax)
    xkeep <- c(ax, cx, x1, x2)
    fkeep <- c(fa, fc, f1, f2)
    if (plot.it == TRUE) {
        plot(c(ax, bx, cx), c(fa, fb, fc))
        text(c(x1, x2), c(f1, f2), lab = c("1", "2"))
    }
    cnt <- 3
 
    while ((abs(x3 - x0) > tol * (abs(x1) + abs(x2))) && its <= maxits) {
        if (verbose > 0) {
            cat("x0=", x0, "x1=", x1, "x2=", x2, "x3=", x3, "\n")
            cat("f1=", f1, "f2=", f2, "\n")
        }
        if (f2 < f1) {
            x0 <- x1
            x1 <- x2
            x2 <- round(R * x1 + C * x3)
            f1 <- f2
            f2.ews <- ewspec3(x = x, filter.number = filter.number, 
                family = family, smooth.dev = smooth.dev, AutoReflect = AutoReflect, 
                WPsmooth.type = "RM", binwidth = x2)$S
            f2 <- specerr(f2.ews, EWS, Jmin:Jmax)
            if (verbose == 2) {
                cat("SSQ: ", signif(f2, 3), "\n")
            }
            else if (verbose == 1) 
                cat(".")
            xkeep <- c(xkeep, x2)
            fkeep <- c(fkeep, f2)
            if (plot.it == TRUE) 
                text(x2, f2, lab = as.character(cnt))
            cnt <- cnt + 1
        }
        else {
            x3 <- x2
            x2 <- x1
            x1 <- round(R * x2 + C * x0)
            f2 <- f1
            f1.ews <- ewspec3(x = x, filter.number = filter.number, 
                family = family, smooth.dev = smooth.dev, AutoReflect = AutoReflect, 
                WPsmooth.type = "RM", binwidth = x1)$S
            f1 <- specerr(f1.ews, EWS, Jmin:Jmax)
            if (verbose == 2) 
                cat("SSQ: ", signif(f1, 3), "\n")
            else if (verbose == 1) 
                cat(".")
            xkeep <- c(xkeep, x1)
            fkeep <- c(fkeep, f1)
            if (plot.it == TRUE) 
                text(x1, f1, lab = as.character(cnt))
            cnt <- cnt + 1
        }
	its <- its + 1
    }
    if (its > maxits)
	warning(paste("AutoBestBW: exceeded maximum number of iterations: ", maxits))


    if (f1 < f2) 
        ans <- x1
    else ans <- x2

    if (ReturnAll == FALSE)
	    return(ans)
    else	{
	l <- list(ans=ans, EWS.wavelet=EWS, EWS.linear=f1.ews)
        return(l)
	}
}
