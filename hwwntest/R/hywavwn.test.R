hywavwn.test <-
function (x, away.from="standard", lowlev = 0, plot.it = FALSE,
	stopeveryscale = FALSE, filter.number=10, family="DaubExPhase",
	mc.method = p.adjust.methods, verbose=FALSE, 
	n.cdf.grid = 1000, mac.spread = 10) 
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

#
# Stuff for both methods
#
    mc.method <- match.arg(mc.method)
    FT = abs(fft(data)/sqrt(N))^2
    periodogram = FT[2:(N/2 + 1)]
    sigsq <- var(x)
    periodogram <- periodogram/sigsq
#
# Stuff for general wavelet for coarser scales
#
    wdecomp.gw = wd(periodogram, filter.number = filter.number,
	family = family, type = "wavelet", bc = "periodic",
	verbose = FALSE, min.scale = 0)
    nlev.gw <- nlevelsWT(wdecomp.gw)
    p.val.collector.gw <- NULL
    for (i in (nlev.gw - 1-away.from):lowlev) {
	if (verbose==TRUE)
		cat("General level: ", i, "\n")
        #m <- 2^(nlev.gw - i - 1)
        wv.gw <- accessD(wdecomp.gw, level = i)
        range.wv.gw <- range(wv.gw)
        fifteen.gw <- (range.wv.gw[2] - range.wv.gw[1]) * 0.3
        range.wv.gw[1] <- range.wv.gw[1] - fifteen.gw
        range.wv.gw[2] <- range.wv.gw[2] + fifteen.gw
        if (plot.it == TRUE) {
            if (length(wv.gw) <= 2) {
                cat("Can't do density plot with less than 3 obs\n")
            }
            else {
                plot(density(wv.gw), main = paste("Wavelet coefficients. Level: ", 
                  i))
                xx <- seq(from = -mac.spread, to = mac.spread, 
                  length = 100)
                yy <- dnorm(xx)
                lines(xx, yy, col = 2)
                if (stopeveryscale == TRUE) 
                  scan()
            }
        }
        if (plot.it == TRUE) {
            plot(ecdf(wv.gw))
            lines(xx, pnorm(xx), col = 2)
            if (stopeveryscale == TRUE) 
                scan()
        }
        wv.cdf.gw <- pnorm(abs(wv.gw))
        wv.pval.gw <- 2 * (1 - wv.cdf.gw)
        if (plot.it == TRUE) {
            hist(wv.pval.gw, main = "Wavelet Coefficient p-values")
            if (stopeveryscale == TRUE) 
                scan()
        }
        p.val.collector.gw <- c(p.val.collector.gw, wv.pval.gw)
    }

    # Don't want to execute following lines from genwwn as we're
    # going to operate on and return combined p-values
    #p.val.adjust.gw <- p.adjust(p.val.collector.gw, method = mc.method)
    #ll <- list(p.val.collector = p.val.collector, p.val.adjust = p.val.adjust, 
    #    p.value = min(p.val.adjust), method = "Wavelet Test of White Noise")
    #class(ll) <- c("htest")
    #return(ll)

#
#   Stuff for Haar wavelet for finest scales.
#

    wdecomp.hw = wd(periodogram, filter.number = 1, family = "DaubExPhase", 
        type = "wavelet", bc = "periodic", verbose = FALSE, min.scale = 0)
    nlev.hw <- nlevelsWT(wdecomp.hw)
    p.val.collector.hw <- NULL
    for (i in (nlev.hw - 1):(nlev.hw - away.from)) {
	if (verbose==TRUE)
		cat("Haar level: ", i, "\n")
        m <- 2^(nlev.hw - i - 1)
        wv.hw <- accessD(wdecomp.hw, level = i)
        range.wv.hw <- range(wv.hw)
        fifteen.hw <- (range.wv.hw[2] - range.wv.hw[1]) * 0.3
        range.wv.hw[1] <- range.wv.hw[1] - fifteen.hw
        range.wv.hw[2] <- range.wv.hw[2] + fifteen.hw
        if (plot.it == TRUE) {
            if (length(wv.hw) <= 2) {
                cat("Can't do density plot with less than 3 obs\n")
            }
            else {
                plot(density(wv.hw), main = paste("Haar Wavelet coefficients. Level: ", 
                  m))
                xx <- seq(from = -mac.spread, to = mac.spread, 
                  length = 100)
                yy <- Macdonald(xx, m = m)
                lines(xx, yy, col = 2)
                if (stopeveryscale == TRUE) 
                  scan()
            }
        }
        xx <- seq(from = -mac.spread, to = mac.spread, length = n.cdf.grid)
        mac.dens.xx <- Macdonald(xx, m = m)
        cdf.mac.hw = approxfun(x = xx, y = cumsum(mac.dens.xx)/sum(mac.dens.xx), 
            yleft = 0, yright = 1)
        if (plot.it == TRUE) {
            plot(ecdf(wv.hw))
            lines(xx, cdf.mac.hw(xx), col = 2)
            if (stopeveryscale == TRUE) 
                scan()
        }
	wv.cdf.hw <- cdf.mac.hw(abs(wv.hw))
        wv.pval.hw <- 2 * (1 - wv.cdf.hw)
        if (plot.it == TRUE) {
            hist(wv.pval.hw, main = "Wavelet Coefficient p-values")
            if (stopeveryscale == TRUE) 
                scan()
        }
        p.val.collector.hw <- c(p.val.collector.hw, wv.pval.hw)
    }

#
# Put general wavelet and Haar wavelet pvalues together
#

    p.val.collector <- c(p.val.collector.gw, p.val.collector.hw)

    p.val.adjust <- p.adjust(p.val.collector, method = mc.method)

    ll <- list(p.val.collector = p.val.collector, p.val.adjust = p.val.adjust, 
        p.value = min(p.val.adjust), method = "Hybrid Wavelet Test of White Noise", p.val.collector.hw=p.val.collector.hw, p.val.collector.gw=p.val.collector.gw)
    class(ll) <- c("htest")
    return(ll)
}
