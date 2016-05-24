genwwn.test <-
function (x, away.from="standard", lowlev = 0, plot.it = FALSE, stopeveryscale = FALSE, 
    filter.number=10, family="DaubExPhase",
    mc.method = p.adjust.methods, mac.spread = 10, verbose=FALSE) 
{
    data = x
    N = length(data)
    if (is.na(J.N <- IsPowerOfTwo(N))) 
        stop("Data has to be of length a power of two")
    if (N < 16)
	stop("Need time series to be of length 16 or greater")
#
# Away.from can be an integer (keep away from the away.from finest scales,
# or the text string "standard" where we internally set the away.from based
# on the length of the data
#
    if (!is.numeric(away.from))	{
	if (away.from != "standard")
		stop("away.from argument has to be positive integer or the string `standard'")
	if (N==16 || N==32)
		away.from <- 2
	else if (N > 32 && N <= 1024)  
		away.from <- J.N - 4
	else if (N > 1024)
		away.from <- J.N - 5

	if (verbose==TRUE)
		cat("N is: ", N, " J.N is: ", J.N, " away.from is: ", away.from, "\n")
	}

    mc.method <- match.arg(mc.method)
    FT = abs(fft(data)/sqrt(N))^2
    periodogram = FT[2:(N/2 + 1)]
    sigsq <- var(x)
    periodogram <- periodogram/sigsq
    wdecomp = wd(periodogram, filter.number = filter.number,
	family = family,
        type = "wavelet", bc = "periodic", verbose = FALSE, min.scale = 0)
    nlev <- nlevelsWT(wdecomp)
    p.val.collector <- NULL
    for (i in (nlev - 1-away.from):lowlev) {
        m <- 2^(nlev - i - 1)
        wv <- accessD(wdecomp, level = i)
        range.wv <- range(wv)
        fifteen <- (range.wv[2] - range.wv[1]) * 0.3
        range.wv[1] <- range.wv[1] - fifteen
        range.wv[2] <- range.wv[2] + fifteen
        if (plot.it == TRUE) {
            if (length(wv) <= 2) {
                cat("Can't do density plot with less than 3 obs\n")
            }
            else {
                plot(density(wv), main = paste("Wavelet coefficients. Level: ", 
                  m))
                xx <- seq(from = -mac.spread, to = mac.spread, 
                  length = 100)
                yy <- dnorm(xx)
                lines(xx, yy, col = 2)
                if (stopeveryscale == TRUE) 
                  scan()
            }
        }
        if (plot.it == TRUE) {
            plot(ecdf(wv))
            lines(xx, pnorm(xx), col = 2)
            if (stopeveryscale == TRUE) 
                scan()
        }
        wv.cdf <- pnorm(abs(wv))
        wv.pval <- 2 * (1 - wv.cdf)
        if (plot.it == TRUE) {
            hist(wv.pval, main = "Wavelet Coefficient p-values")
            if (stopeveryscale == TRUE) 
                scan()
        }
        p.val.collector <- c(p.val.collector, wv.pval)
    }
    p.val.adjust <- p.adjust(p.val.collector, method = mc.method)
    ll <- list(p.val.collector = p.val.collector, p.val.adjust = p.val.adjust, 
        p.value = min(p.val.adjust), method = "Wavelet Test of White Noise")
    class(ll) <- c("htest")
    return(ll)
}
