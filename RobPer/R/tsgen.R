tsgen <- function(ttype, ytype, pf, redpart, s.outlier.fraction=0, interval, npoints, ncycles, ps, SNR, alpha=1.5) {
    # alpha   : alpha for red noise
    # interval: logical: whether the y-observations in an interval of length 3ps should be replaced by a peak or not
    # npoints : number of observations
    # ncycles: length of observation interval, in cycles of length ps
    # pf     : fluctuation period
    # ps     : sampling period
    # redpart           : fraction in [0,1[ for red noise of the overall noise dispersion
    # s.outlier.fraction: A fraction of measurement accuracies that is replaced by outliers	
    # SNR               : Signal- to- noise ratio: Var(sig)/var(noise)
    # ttype: see sampler.R
    # ytype: see signalgen.R
  
    # Sample time points
    tt <- sampler(ttype, npoints, ncycles, ps)
	
    # calculate signal
    sig <- signalgen(tt, ytype, pf)
	
    # Add noise and scale signal to the right SNR
    temp <- lc_noise(tt,sig, SNR, redpart, alpha=alpha)
    y <- temp$y
    s <- temp$s
	
    # replace s-observations by tiny outliers or include a peak
    temp <- disturber(tt,y,s,ps, s.outlier.fraction=s.outlier.fraction, interval)
    y <- temp$y
    s <- temp$s
	
    return(cbind(tt, y, s))
}
