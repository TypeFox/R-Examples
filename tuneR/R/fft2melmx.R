# This code is based on the Matlab implementations of PLP and Rasta
# feature calculations by Daniel P. W. Ellis of Columbia University /
# International Computer Science Institute.  For more details, see:
# http://www.ee.columbia.edu/~dpwe/resources/matlab/rastamat/

fft2melmx <- function(nfft, sr=8000, nfilts=40, width=1.0, minfreq=0, maxfreq=sr/2,
                htkmel=FALSE, constamp=FALSE){

    if(!(sr==as.integer(sr) && nfft==as.integer(nfft)) || sr <= 0 || nfft <= 0)
      stop("'sr' and 'nfft' have to be positive integers")

    if(!is.null(nfilts) && !(nfilts==as.integer(nfilts) && nfilts > 0))
      stop("'nfilts' has to be positive and integer valued")

    fftfreqs <- (0:(nfft-1))/nfft * sr

    minmel <- hz2mel(f=minfreq, htk=htkmel)
    maxmel <- hz2mel(f=maxfreq, htk=htkmel)
    # Frequency of each FFT bin in Mel
    binfreqs <- mel2hz(z=minmel + (0:(nfilts+1))/(nfilts+1) * (maxmel-minmel), htk=htkmel)

    binbin <- round(binfreqs/sr * (nfft-1))

    wtscalc <- function(i, binfreqs=binfreqs){
        fs <- binfreqs[i + c(0, 1, 2)]
        # Scale by width
        fs <- fs[2] + width * (fs - fs[2])

        # Calculate slopes
        loslope <- (fftfreqs - fs[1])/(fs[2] - fs[1])
        hislope <- (fs[3] - fftfreqs)/(fs[3] - fs[2])

        return(pmax(0, pmin(loslope, hislope)))
    }
    wts <- t(sapply(seq(nfilts), function(x) wtscalc(i=x, binfreqs=binfreqs)))

    if(!constamp){
        # Scale to be approx constant E (Slaney-style mel)
        wts <- diag(2/(binfreqs[2+(1:nfilts)] - binfreqs[1:nfilts])) %*% wts
    }

    # Ensure 2nd half of FFT ist zero
#    wts[,(nfft/2 + 1):nfft] <- 0

    return(list(wts=wts, binfreqs=binfreqs))
}
