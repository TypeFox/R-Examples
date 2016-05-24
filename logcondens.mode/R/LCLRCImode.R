## Charles Doss
## function which computes confidence interval for mode based on
## the log-concave likelihood ratio.


#### Now use sysdata.rda
##data(TLLRdistn)
#### Now just list 'distr' in Imports in DESCRIPTION, and access via "::".
##library(distr)
## this is just in sysdata.rda now
##LCTLLRdistn <- distr::DiscreteDistribution(supp=TLLRsUnq$TLLRs, prob=TLLRsUnq$ws / sum(TLLRsUnq$ws))

LRmodeTest <-  
  function(mode, x, xgrid=NULL, w=NA,
           nn=length(x),
           alpha,prec=1e-10,print=FALSE){
    ## redundant, because this check is done in activeSetLogCon and in
    ## activeSetLogCon.mode
    if ((!identical(xgrid, NULL) & (!identical(w, NA)))) {
      stop("If w != NA then xgrid must be NULL!\n")
    }
    if ((!is.numeric(w) || !is.null(xgrid)) && (nn != length(x))) {
      stop("If nn!=length(x), then w must be given a value and xgrid must be NULL.")
    }
    res.UC <- activeSetLogCon(x=x,w=w,
                              xgrid=xgrid,
                              prec=prec,print=print) ##UnConstrained
    res.MC <- activeSetLogCon.mode(x=x,
                                   xgrid=xgrid,
                                   w=w,
                                   mode=mode,
                                   prec =prec,
                                   print=print) ##Mode Constrained
    LL <- 2 * nn * (res.UC$L - res.MC$L); 
    ##CC <- logcondens.mode:::LCTLLRdistn@q(1-alpha)
    CC <- LCTLLRdistn@q(1-alpha)
    if (LL <= CC) TRUE
    else FALSE
  }


## Note that if activeSetLogCon and activeSetLogCon.mode took starting values
## for the density estimates, this function could be sped up considerably.

## And, of course, there is no need to re-compute res.UC each iteration. And
## I then compute res.UC _again_ to get the mode estimate.  So this is
## already inefficient by a factor of 3.

LCLRCImode <- 
  function(x,
           xgrid=NULL,
           w=NA,
           nn=length(x), ##nn is used if w is not NA and xgrid is NULL
           alpha=.05,
           prec=1e-10,
           CIprec=1e-4,
           print=F){


    ##print("TEST THIS; should be named 'lclrcimode.weights") 

    ## redundant, this check is done by LRmodeTest, but.
  if ((!identical(xgrid, NULL) & (!identical(w, NA)))) {
    stop("If w != NA then xgrid must be NULL!\n")
  }
  if ((!is.numeric(w) || !is.null(xgrid)) && (nn != length(x))) {
    stop("If nn!=length(x), then w must be given a value and xgrid must be NULL.")
  }
  if (!is.numeric(nn) || length(nn) > 1){
    stop(paste("nn was passed a value of", nn,
               ", which is not acceptable; should be a numeric of length 1",
               sep=""))
  }

  if (!is.numeric(w)){  ##xgrid may or may not be NULL, preProcess handles either 
    tmp <-  preProcess(x=x, xgrid=xgrid)
    x <- tmp$x
    w <- tmp$w
    nn <- tmp$n ##the length of the original x
  }
  myLRmodeTest <- function(mm){LRmodeTest(mode=mm,
                                          x=x,
                                          ## xgrid, ##preprocess
                                          w=w,
                                          nn=nn,
                                          alpha=alpha,prec=prec,print=print)}
  MLE.UC <- activeSetLogCon(x=x,
                            w=w,
                            ##xgrid=xgrid,
                            prec=prec,
                            print=print) ##UnConstrained
  mhat <- MLE.UC$dlcMode$val
  if (!myLRmodeTest(mhat)) return(numeric(0)) ## then something weird happening. 

  {## now get outer limits of the interval
    supp.len <- x[length(x)]-x[1];
    L.start <- x[1] ## L.start, R.start are starting values for left,right endpoints of interval.
    R.start <- x[length(x)] ##they usually are x[1], x[length(x)].
    ## arbitrary: linearly decrease length until exclude endpoint
    while (myLRmodeTest(L.start))
      L.start <- L.start - supp.len 
    while (myLRmodeTest(R.start))
      R.start <- R.start + supp.len
  }
  ## now use binary search starting from L.start, R.start
  ## note that test rejects L.start and R.start
  {
    L <- L.start
    R <- mhat ## is always right of R.start
    while(R-L > CIprec){ ##R is in, L is out
      mid <- (R+L)/2
      if (myLRmodeTest(mid)) R <- mid
      else L <- mid
    }
  }
  Lpt <- L
  ## Do right side
  ##  if (myLRmodeTest(x[length(x)])) R <- x[length(x)]   else
  {
    L <- mhat
    R <- R.start
    while (R-L > CIprec){ ## R is out L is in
      mid <- (R+L)/2
      if (myLRmodeTest(mid)) L <- mid
      else R <- mid
    }
  }
  Rpt <- R
  return(c(Lpt,Rpt))
}



