## $Id: baseline.rfbaseline.R 182 2011-01-09 21:05:18Z kristl $
baseline.rfbaseline <- function(spectra, span=2/3, NoXP=NULL, maxit=c(2,2), b=3.5,
                                weight=NULL, Scale=function(r) median(abs(r))/0.6745,
                                delta=NULL, SORT=FALSE, DOT=FALSE, init=NULL){
  
  if(requireNamespace("IDPmisc", quietly = TRUE)){
    
    np <- dim(spectra)
    baseline  <- matrix(0,np[1],np[2])
    X <- 1:np[2]
    
    # Sends one spectrum at the time to rfbaseline
    for(i in 1:np[1]){
      rbe <- IDPmisc::rfbaseline(x=X, y=spectra[i,], span=span, NoXP=NoXP, maxit=maxit, b=b,
                        weight=weight, Scale=Scale,
                        delta=delta, SORT=SORT, DOT=DOT, init=init)
      baseline[i,]  <- rbe$fit
    }
    return(list(baseline = baseline, corrected = spectra - baseline))
  } else {
    warning('Package IDPmisc not installed')
    return(list())
  }
}
