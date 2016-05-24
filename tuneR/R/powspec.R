# This code is based on the Matlab implementations of PLP and Rasta
# feature calculations by Daniel P. W. Ellis of Columbia University /
# International Computer Science Institute.  For more details, see:
# http://www.ee.columbia.edu/~dpwe/resources/matlab/rastamat/

powspec <- function(x, sr=8000, wintime=0.025, steptime=0.010, dither=FALSE){
    if((!is.numeric(x)) || (!is.null(dim(x))))
      stop("'x' has to be a numeric vector.")

    if(!(sr==as.integer(sr) && is.numeric(wintime) && is.numeric(steptime)) ||
      sr <= 0 || wintime <= 0 || steptime <= 0)
      stop("'sr', 'wintime' and 'steptime' have to be positive; 'sr' also
        integer valued")

    winpts <- round(wintime * sr)
    steppts <- round(steptime * sr)

    nfft <- 2^(ceiling(log(winpts)/log(2)))
    window <- hamming(winpts)
    noverlap <- winpts - steppts

    y <- abs(specgram(x=x, n=nfft, Fs=sr, window=window, overlap=noverlap)$S)^2

    # Avoid digital zero
    if(dither){
      y <- y + winpts
    }

    return(y)
}
