lc_noise <- function(tt,sig, SNR, redpart, alpha=1.5) {
    if(redpart>=1)              stop("redpart has to lie in [0,1[")
    if(length(tt)!=length(sig)) stop("vectors tt and sig must have the same length")
	
    constant.signal<-var(sig)==0
    npoints  <- length(tt)
    redtimes <- (redpart)/(1-redpart) # redtimes*var(white)= var(red)
    # Measurement accuracies
    s <- rgamma(npoints, shape=3, rate=10)   # mean of distribution: shape/rate
  
    # noise
    white <- rnorm(npoints, sd=s)
    if(redpart>0) {
        red_ <- TK95_uneq(tt, alpha)
        red <- red_*sd(white)/sd(red_)*sqrt(redtimes)
    } else {
        red<- 0
    }
    noise <- white+ red
	
    # Scale the noise (new in order to get the right SNR if the signal varies at all
    if (constant.signal) c <- 1 else c <- sqrt(var(sig)/ (SNR*var(noise)) )
    noise <- c*noise
	
    y<- sig+noise	
    return(list(y=y, s=s))
}
