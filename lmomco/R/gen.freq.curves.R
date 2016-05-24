"gen.freq.curves" <-
function(n, para, F=NULL, nsim=10, callplot=TRUE, aslog=FALSE, asprob=FALSE,
         showsample=FALSE, showparent=FALSE, lowerCI=NA, upperCI=NA, FCI=NA,...) {

  if(! are.par.valid(para)) return()
  type  <- para$type

  if(is.null(F)) F <- nonexceeds()
  if(! check.fs(F)) return()

  if(! is.na(upperCI) && length(upperCI) > 1) {
     warning("upperCI argument is a vector, using only first value")
     upperCI <- upperCI[1]
  }
  if(! is.na(lowerCI) && length(lowerCI) > 1) {
     warning("lowerCI argument is a vector, using only first value")
     lowerCI <- lowerCI[1]
  }
  if(! is.na(FCI) && length(FCI)) {
     warning("FCI argument is a vector, using only first value")
     lowerCI <- lowerCI[1]
  }

  plotF  <- F
  xlabel <- 'NONEXCEEDANCE PROBABILITY'
  ylabel <- 'QUANTILE'
  if(asprob == TRUE) {
    plotF  <- qnorm(F)
    xlabel <- 'STANDARD NORMAL QUANTILE'
  }

  if(callplot == TRUE) {
    Q <- par2qua(F,para)
    if(aslog == TRUE) {
      Q <- log10(Q)
      ylabel <- 'log10(QUANTILE)'
    }
    plot(plotF,Q,type='n',xlab=xlabel,ylab=ylabel, ...)
  }

  count.lowerCI <- ifelse(is.na(lowerCI), NA, 0)
  count.upperCI <- ifelse(is.na(upperCI), NA, 0)

  count <- 0
  while(count < nsim) {
    sX   <- rlmomco(n,para)
    sLMR <- lmoms(sX)
    if(! are.lmom.valid(sLMR)) {
      # yet another ad hoc solution for bizarre simulated values
      next
    }
    sPAR <- lmom2par(sLMR,type=type)
    if(is.null(sPAR) || ( length(sPAR$ifail) == 1
                                   &&
                            sPAR$ifail  != 0 ) ) {
      # The ifail is suitable for kappa (kap) and wakeby (wak)
      # distributions.
      # if the parameters could not be solved
      # for the desired distribution of the
      # parent---just do it again---ad hoc.
      next
    }
    count  <- count + 1

    QCI <- ifelse(is.na(FCI), NA, par2qua(FCI, sPAR))
    if(! is.na(QCI)) {
       if(!is.na(upperCI) && QCI >= upperCI) count.upperCI <- count.upperCI + 1
       if(!is.na(lowerCI) && QCI <= lowerCI) count.lowerCI <- count.lowerCI + 1
    }
    freqcurve <- par2qua(F,sPAR)
    if(aslog == TRUE) freqcurve <- log10(freqcurve)
    lines(plotF,freqcurve,...)
    if(showsample == TRUE) {
      sX <- sort(sX)
      plotting.position <- pp(sX)
      if(asprob == TRUE) plotting.position <- qnorm(plotting.position)
      if(aslog == TRUE)  sX <- log10(sX)
      points(plotting.position,sX,...)
    }
  }
  if(showparent == TRUE) {
    parent <- par2qua(F,para)
    if(aslog == TRUE) parent <- log10(parent)
    lines(plotF,parent,lwd=3)
  }
  z <- list(nonexceedance.probability = FCI,
            count.above.upperCI       = count.upperCI,
            count.below.lowerCI       = count.lowerCI,
            count.valid.simulations   = count)
  return(z)
}
