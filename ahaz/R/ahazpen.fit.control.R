"ahazpen.fit.control"<-function(thresh=1e-5,maxit=100000,...)
  {
    ## Purpose: (Internal) control function for CCD algorithm;
    ##          to be used only in ahazpen calls
    ## ----------------------------------------------------------------------
    ## Arguments:
    ##   thresh         : relative threshold for declaring convergence
    ##   maxit          : maximal number of iterations - may be necessary to
    ##                    increase when close-to-saturated solution required
    ## ----------------------------------------------------------------------
    ## Author: Anders Gorst-Rasmussen
    
    if(!is.numeric(thresh) || thresh<=0)
      stop("value of 'thresh' must be > 0")
    if(!is.numeric(maxit) || maxit <=0)
      stop("maximum number of iterations must be >0")
    list(thresh=thresh, maxit=maxit)
  }
