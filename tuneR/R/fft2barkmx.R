# This code is based on the Matlab implementations of PLP and Rasta
# feature calculations by Daniel P. W. Ellis of Columbia University /
# International Computer Science Institute.  For more details, see:
# http://www.ee.columbia.edu/~dpwe/resources/matlab/rastamat/

fft2barkmx <- function(nfft, sr=8000, nfilts=NULL, width=1.0, minfreq=0, maxfreq=sr/2){

    if(!(sr==as.integer(sr) && nfft==as.integer(nfft)) || sr <= 0 || nfft <= 0)
      stop("'sr' and 'nfft' have to be positive integers")

    if(!is.null(nfilts) && !(nfilts==as.integer(nfilts) && nfilts > 0))
      stop("'nfilts' has to be positive and integer valued")

    min_bark <- hz2bark(f=minfreq)
    nyqbark <- hz2bark(f=maxfreq) - min_bark

    if(is.null(nfilts)){
        nfilts <- ceiling(nyqbark) + 1
    }

    wts <- matrix(0, nrow=nfilts, ncol=nfft)

    # Filterspacing in Bark
    step_barks <- nyqbark/(nfilts-1)

    # Frequency of each FFT bin in Bark
    binbarks <- hz2bark( f=(0:(nfft/2)) * sr/nfft )

    wtscalc <- function(i, min_bark=min_bark, step_barks=step_barks,
                binbarks=binbarks){
        f_bark_mid <- min_bark + (i-1) * step_barks
        # Linear slopes in logarithmic space
        lof <- (binbarks - f_bark_mid)/width - 0.5
        hif <- (binbarks - f_bark_mid)/width + 0.5
        return(10^(pmin(0, pmin(hif, -2.5 * lof))))
    }
    wts[,1:(nfft/2+1)] <- t(sapply(seq(nfilts), function(x) wtscalc(i=x,
                min_bark=min_bark, step_barks=step_barks, binbarks=binbarks)))
    return(list(wts=wts))
}
