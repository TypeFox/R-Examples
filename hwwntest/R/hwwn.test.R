hwwn.test <-
function (x, lowlev = 0, plot.it = FALSE, stopeveryscale = FALSE, 
    n.cdf.grid = 1000, mc.method = p.adjust.methods, mac.spread = 10) 
{
    mc.method <- match.arg(mc.method)
    data = x
    N = length(data)
    if (is.na(IsPowerOfTwo(N))) 
        stop("Data has to be of length a power of two")
    FT = abs(fft(data)/sqrt(N))^2
    periodogram = FT[2:(N/2 + 1)]
    sigsq <- var(x)
    periodogram <- periodogram/sigsq
    wdecomp = wd(periodogram, filter.number = 1, family = "DaubExPhase", 
        type = "wavelet", bc = "periodic", verbose = FALSE, min.scale = 0)
    nlev <- nlevelsWT(wdecomp)
    p.val.collector <- NULL
    for (i in (nlev - 1):lowlev) {
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
                yy <- Macdonald(xx, m = m)
                lines(xx, yy, col = 2)
                if (stopeveryscale == TRUE) 
                  scan()
            }
        }
        xx <- seq(from = -mac.spread, to = mac.spread, length = n.cdf.grid)
        mac.dens.xx <- Macdonald(xx, m = m)
        cdf.mac = approxfun(x = xx, y = cumsum(mac.dens.xx)/sum(mac.dens.xx), 
            yleft = 0, yright = 1)
        if (plot.it == TRUE) {
            plot(ecdf(wv))
            lines(xx, cdf.mac(xx), col = 2)
            if (stopeveryscale == TRUE) 
                scan()
        }
        wv.cdf <- cdf.mac(abs(wv))
        wv.pval <- 2*(1-wv.cdf)
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
