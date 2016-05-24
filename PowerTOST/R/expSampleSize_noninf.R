#------------------------------------------------------------------------------
# Sample size based on 'expected' power
# taking into account the uncertainty of CV for the non-inferiority
# t-test
# 
# Author: dlabes
#------------------------------------------------------------------------------
# sample size start for Julious 'expected' power
.expsampleN0.noninf <- function(alpha, targetpower, lmargin, diffm, 
                                se, dfse, steps=2, bk=2)
{
  Z1 <- qnorm(1-alpha)
  tinv <- qt(targetpower, dfse, Z1)
  # factor 2 in Julious = bk
  n0 <- ceiling(bk*(se*tinv/(diffm - lmargin))^2)
  #make an even multiple of step (=2 in case of 2x2 cross-over)
  n0 <- steps*trunc(n0/steps)
  #if (n0<4) n0 <- 4   # minimum sample size will be tested outside
  
  return(n0)
}
#------------------------------------------------------------------------------
# Sample size for a desired "expected" power according to Julious: 
# see known.designs() for covered experimental designs
# Only for log-transformed data
# CV and dfCV can be vectors, if then a pooled CV, df will be calculated
expsampleN.noninf <- function(alpha=0.025, targetpower=0.8, logscale=TRUE,
                              theta0, margin, CV, dfCV, design="2x2", robust=FALSE, 
                              method=c("exact", "approx"), print=TRUE, 
                              details=FALSE, imax=100)
{
  #number of the design and check if design is implemented
  d.no <- .design.no(design)
  if (is.na(d.no)) stop("Design ",design," not known!", call.=FALSE)
  
  # design characteristics
  ades   <- .design.props(d.no)
  d.name <- ades$name  # nice name of design
  # get the df for the design as an unevaluated expression
  dfe    <- .design.df(ades, robust=robust)
  steps  <- ades$steps	# stepsize for sample size search
  bk     <- ades$bk     # get design constant
  # minimum sample size
  df <- 0; n <- 0
  while (df<1){
    n  <- n + 1
    df <- eval(dfe)
  }
  # make a multiple of steps
  nmin <- as.integer(steps*trunc(n/steps)) 
  nmin <- nmin + steps*(nmin<n)
  
  if (missing(CV) | missing(dfCV)) {
    stop("CV and its df must be given!", call.=FALSE)
  }
  # calculate pooled data if CV and dfCV are vectors
  if (length(CV)>1){
    if (length(dfCV)!=length(CV)) {
      stop("CV and df must have equal number of entries!", call.=FALSE)
    }
    dfse <- sum(dfCV)
    if (logscale) {
      CVp  <- CV2se(CV)^2 #need s-squared
    } else {
      CVp <- CV^2
    }  
    CVp  <- CVp * dfCV
    CVp  <- sum(CVp)/dfse
    CVp  <- sqrt(CVp)
    if(logscale) CVp  <- se2CV(CVp) 
  } else {
    dfse <- dfCV
    CVp  <- CV
  }
  
  # check method
  method <- tolower(method)
  method <- match.arg(method)
  
  # print the configuration:
  if (print) {
    cat("\n+++++++++++ Non-inferiority test +++++++++++\n")
    cat("     Sample size est. with uncertain CV\n")
    cat("--------------------------------------------\n")
    cat("Study design: ",d.name,"\n")
    if (details) { 
      cat("Design characteristics:\n")
      if (robust & (ades$df2 != ades$df)) {
        cat("df = ",ades$df2," (robust)", sep="") 
      } else cat("df = ",ades$df, sep="")
      cat(", design const. = ", bk, ", step = ", steps,"\n\n",sep="")
    }     
  }

  # handle the log transformation
  if (logscale) {
    if (missing(margin)) margin <- 0.8
    if (missing(theta0)) theta0 <- 0.95
    if ( (theta0<=margin) & (margin<1) ) {
      stop("Null ratio ",theta0," must be above margin ",margin,"!", 
          call.=FALSE)
    }
    if ( (theta0>=margin) & (margin>1) ) {
      stop("Null ratio ",theta0," must be below margin ",margin,"!", 
          call.=FALSE)
    }
    lmargin <- log(margin)
    diffm   <- log(theta0)
    se      <- CV2se(CVp)
    if (print) cat("log-transformed data (multiplicative model)\n\n")
  } else {
    if (missing(margin)) margin <- -0.2
    if (missing(theta0)) theta0 <- -0.05
    if ( (theta0<=margin) & (margin<0) ) {
      stop("Null diff. ",theta0," must be above margin ",margin,"!", call.=FALSE)
    }
    if ( (theta0>=margin) & (margin>0) ) {
      stop("Null diff. ",theta0," must be below margin ",margin,"!", call.=FALSE)
    }
    lmargin <- margin
    diffm   <- theta0
    se      <- CVp
    if (print) cat("untransformed data (additive model)\n\n")
  }
  
  if (print) {
    cat("alpha = ",alpha,", target power = ", targetpower,"\n", sep="")
    cat("Non-inf. margin    = ", margin, "\n", sep="")
    cat("Null (true) ratio  = ",theta0,"\n", sep="")
    if (length(CV)>1){
      cat("Variability data\n")
      print(data.frame(CV=CV,df=dfCV), row.names = FALSE)
      cat("CV(pooled)         = ", CVp, " with ", dfse," df\n", sep="")
    } else {
      cat("CV                 = ", CVp, " with ", dfse," df\n", sep="")
    }   
  }
  
  #start value from large sample approx. 
  n   <- .expsampleN0.noninf(alpha, targetpower, lmargin, diffm, 
                             se, dfse, steps, bk)
  if(n<nmin) n <- nmin
  df  <- eval(dfe)
  #browser()
  pow <- .exppower.noninf(alpha=alpha, lmargin=lmargin, ldiff=diffm, 
                          se.fac=sqrt(bk/n), se=se, dfCV=dfse, df=df, 
                          method=method) 
  if (details) {
    cat("\nSample size search (ntotal)\n")
    cat(" n    exp. power\n")
    # do not print first too high
    if (pow<=targetpower) cat( n," ", formatC(pow, digits=6, format="f"),"\n")
  }
  # --- loop until power >= target power
  iter <- 0
  # iter>50 is emergency brake
  # this is eventually not necessary, depends on quality of sampleN0
  # in experimentation I have seen max of six steps
  # starting with too high power should be rare sinsce the large sample
  # approximation should give too low sample size
  while (pow>targetpower) {
    if (n<=nmin) { # min number
      if (print & iter==0) cat( n," ", formatC(pow, digits=6, format="f"),"\n")
      break
    }
    n    <- n-steps     # step down 
    iter <- iter+1
    df   <- eval(dfe)
    pow <- .exppower.noninf(alpha=alpha, lmargin=lmargin, ldiff=diffm, 
                            se.fac=sqrt(bk/n), se=se, dfCV=dfse, df=df, 
                            method=method) 
    # do not print first step down
    if (details) cat( n," ", formatC(pow, digits=6),"\n")
    if (iter>imax) break  
  }
  while (pow<targetpower) {
    n    <- n+steps
    iter <- iter+1
    df   <- eval(dfe)
    pow <- .exppower.noninf(alpha=alpha, lmargin=lmargin, ldiff=diffm, 
                            se.fac=sqrt(bk/n), se=se, dfCV=dfse, df=df, 
                            method=method) 
    if (details) cat( n," ", formatC(pow, digits=6, format="f"),"\n")
    if (iter>imax) break 
  }
  if (pow<targetpower) {
    n <- NA
    if (details) cat("Sample size search failed!\n")
  }
  
  if (print && !details) {
    cat("\nSample size (ntotal)\n")
    # parallel group design is now also handled in terms of ntotal
    #if (d.no == 0) cat("(n is sample size per group)\n") 
    cat(" n    exp. power\n")
    cat( n," ", formatC(pow, digits=6, format="f"),"\n")
    if (is.na(n)) cat("Sample size search failed!\n")
  }
  if (details && print) {
    if (method=="exact") 
      cat("\nExact expected power calculation.\n")
  }
  # always print if approx.
  if (print & (method!="exact")){
    approx <- "Approximate expected power calculation \nacc. to Julious/Owen."
    cat("\n",approx,"\n",sep="")
  } 
  
  if (print) cat("\n")
  
  #return results as data.frame
  res <- data.frame(design=design, alpha=alpha, CV=CV, dfCV=dfse, theta0=theta0, 
                    margin=margin, n=n, power=pow, targetpower=targetpower)
  names(res) <-c("Design","alpha","CV","df of CV","theta0","Margin",
                 "Sample size", "Achieved power", "Target power")
  
  if (print) return(invisible(res)) 
    else return(res)
  
} # end function


