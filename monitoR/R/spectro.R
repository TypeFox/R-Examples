# Based on seewave's spectro function, very few modifications
# Modified: 2012 DEC 22

spectro <-
function (
   wave, 
   wl = 512, 
   wn = "hanning", 
   zp = 0, 
   ovlp = 0, 
   fftw = FALSE, 
   dB = "max0", 
   dBref = NULL) 
{
    if (!is.null(dB) && all(dB != c("max0", "A", "B", "C", "D"))) 
        stop("'dB' has to be one of the following character strings: 'max0', 'A', 'B', 'C' or 'D'")
    f <- wave@samp.rate
    wave <- as.matrix(wave@left)
    n <- nrow(wave)
    #step <- seq(1, n - wl, wl - (ovlp*wl/100)) # This is original, and will just leave out any partial stuff at the end
    step <- seq(1, n - wl + 1, wl - (ovlp*wl/100)) # New version squeezes out one more time bin in some cases
    # Two lines below drop parts of wave that go beyond step vector. I added them, but are they needed? I think so, at least for n, so the X is correct
    # New calculation of n is a bit tricky because of tails with no overlap at beginning and end
    n <- length(step)*wl  - (length(step)-1)*ovlp*wl/100 
    wave <- wave[1:n, ,drop=FALSE]
    z <- stft(wave = wave, f = f, wl = wl, zp = zp, step = step, wn = wn, fftw = fftw)
    X <- seq(0, (n - wl)/f, length.out = length(step)) # X is time, and here is for the left edge of time bins (could be right too, or center even). Note that the left edge of the last time bin is not affected by ovlp.
    Y <- seq((f/1000)/(wl + zp), f/2000, length.out = nrow(z))

    if (!is.null(dB)) {
        if (is.null(dBref)) 
            z <- 20 * log10(z)
        else z <- 20 * log10(z/dBref)
        if (dB == "max0") 
            z <- z
        if (dB == "A") 
            z <- dBweight(Y * 1000, dBref = z)$A
        if (dB == "B") 
            z <- dBweight(Y * 1000, dBref = z)$B
        if (dB == "C") 
            z <- dBweight(Y * 1000, dBref = z)$C
        if (dB == "D") 
            z <- dBweight(Y * 1000, dBref = z)$D
    }

    return(list(time = X, freq = Y, amp = z))
}
