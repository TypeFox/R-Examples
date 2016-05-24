lacf <-
function (x, filter.number = 10, family = "DaubLeAsymm", smooth.dev = var, 
    AutoReflect = TRUE, lag.max = NULL, WPsmooth.type = "RM", 
    binwidth, tol=0.1, maxits=5, ABBverbose=0, verbose=FALSE, ...) 
{
    dsname = deparse(substitute(x))

    if (WPsmooth.type=="RM")	{
	if (missing(binwidth) || binwidth==0)
		binwidth <- AutoBestBW(x=x, filter.number = filter.number,
		family = family, smooth.dev = smooth.dev, 
		AutoReflect = AutoReflect, tol = tol, maxits = maxits, 
		plot.it = FALSE, verbose = ABBverbose)
	if (verbose==TRUE)
		cat("Linear Smoothing. Bandwidth is: ", binwidth, "\n")
	}


    EWS <- ewspec3(x = x, filter.number = filter.number, family = family, 
        smooth.dev = smooth.dev, AutoReflect = AutoReflect, WPsmooth.type = WPsmooth.type, 
        binwidth = binwidth, ...)
    S <- EWS$S
    SmoothWP <- EWS$SmoothWavPer
    J <- S$nlevels
    Smat <- matrix(S$D, nrow = length(x), ncol = J)
    Psi <- PsiJmat(-J, filter.number = filter.number, family = family)
    nc <- ncol(Psi)
    L <- (nc - 1)/2
    dimnames(Psi) <- list(NULL, c(-L:0, 1:L))
    if (is.null(lag.max)) 
        # Old way the.lacf <- Smat %*% Psi[, (L + 1):ncol(Psi)]
	lag.max <- floor(10 * (log10(length(x))))

    if (L + 1 + lag.max > ncol(Psi)) {
        warning(paste("lag.max too high. Have reset it to ", 
                ncol(Psi) - L - 1, ". Higher lags are zero"))
            lag.max <- ncol(Psi) - L - 1
        }
    the.lacf <- Smat %*% Psi[, (L + 1):(L + 1 + lag.max)]
    the.lacor <- sweep(the.lacf, 1, the.lacf[, 1], FUN = "/")
    l <- list(lacf = the.lacf, lacr = the.lacor, name = dsname, 
        date = date(), SmoothWP = SmoothWP, S = S, J = J)
    class(l) <- "lacf"
    return(l)
}
