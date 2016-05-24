# This code is based on the Matlab implementations of PLP and Rasta
# feature calculations by Daniel P. W. Ellis of Columbia University /
# International Computer Science Institute.  For more details, see:
# http://www.ee.columbia.edu/~dpwe/resources/matlab/rastamat/

audspec <- function(pspectrum, sr=16000, nfilts=ceiling(hz2bark(sr/2))+1,
                fbtype=c("bark", "mel", "htkmel", "fcmel"), minfreq=0,
                maxfreq=sr/2, sumpower=TRUE, bwidth=1.0){

    if(!(is.numeric(pspectrum) && is.matrix(pspectrum)))
      stop("'pspectrum' has to be a numeric matrix")

    if(!(sr==as.integer(sr) && sr > 0))
      stop("'sr' has to be a positive integer")

    nfreqs <- nrow(pspectrum)

    nfft <- (nfreqs-1) * 2

    # Construct weight matrix
    fbtype <- match.arg(fbtype)
    wts <- switch(fbtype,
            bark = fft2barkmx(nfft=nfft, sr=sr, nfilts=nfilts, width=bwidth,
                minfreq=minfreq, maxfreq=maxfreq)$wts,
            mel = fft2melmx(nfft=nfft, sr=sr, nfilts=nfilts, width=bwidth,
                minfreq=minfreq, maxfreq=maxfreq)$wts,
            htkmel = fft2melmx(nfft=nfft, sr=sr, nfilts=nfilts, width=bwidth,
                minfreq=minfreq, maxfreq=maxfreq, htkmel=TRUE, constamp=TRUE)$wts,
            fcmel = fft2melmx(nfft=nfft, sr=sr, nfilts=nfilts, width=bwidth,
                minfreq=minfreq, maxfreq=maxfreq, htkmel=TRUE, constamp=FALSE)$wts
    )

    wts <- wts[,1:nfreqs,drop=FALSE]

    # Integrate FFT bins into Mel bins, in abs or abs^2 domain
    if(sumpower){
        aspectrum <- wts %*% pspectrum
    } else {
        aspectrum <- (wts %*% sqrt(pspectrum))^2
    }

    return(list(aspectrum=aspectrum, wts=wts))
}

