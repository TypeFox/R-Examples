#-----------------------------------------------------------------------------
# Sample size calculation based on non-inferiority test
# 
# Author: dlabes
#-----------------------------------------------------------------------------

# start value for sample size search
.sampleN0.noninf <- function(alpha=0.025, targetpower=0.8, lmargin, d0, se, 
                             steps=2, bk=2)
{ 
  n0 <- bk*se^2*(qnorm(targetpower)+ qnorm(1-alpha))^2 / (d0 - lmargin)^2
  n0 <- steps*trunc(n0/steps)
  #if (n0<4) n0 <- 4   # minimum sample size will be tested outside
  return(n0)
}
# --------------------------------------------------------------------------
# Sample size estimation for non-inferiority t-test
sampleN.noninf <- function(alpha=0.025, targetpower=0.8, logscale=TRUE, 
                           margin, theta0, CV, design="2x2", robust=FALSE,
                           details=FALSE, print=TRUE, imax=100)
{ 
  if (missing(CV)) stop("CV must be given!", call.=FALSE)
    
  #number of the design and check
  d.no <- .design.no(design)
  if (is.na(d.no)) stop("Design ",design, " unknown!", call.=FALSE)
  
  # design characteristics
  ades   <- .design.props(d.no)
  d.name <- ades$name   # 'nice' name of design
  # get the df for the design as an unevaluated expression (now with n as var)
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
  if (print) {
    cat("\n++++++++++++ Non-inferiority test +++++++++++++\n")
    cat("            Sample size estimation\n")
    cat("-----------------------------------------------\n")
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
    if(margin<=0) stop("With logscale=TRUE margin must be ratio >0")
    if(theta0<=0) stop("With logscale=TRUE theta0 must be ratio >0")
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
    se      <- CV2se(CV)
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
    se      <- CV
    if (print) cat("untransformed data (additive model)\n\n")
  }
  if (print) {
    cat("alpha = ",alpha,", target power = ", targetpower,"\n", sep="")
    cat("Non-inf. margin   = ", margin, "\n", sep="")
    if (logscale) cat("Null (true) ratio = ",theta0,",  CV = ",CV,"\n", sep="")
    else          cat("Null (true) diff. = ",theta0,",  CV = ",CV,"\n", sep="")
  }
  # start value of 'brute force'
  n   <- .sampleN0.noninf(alpha, targetpower, lmargin, d0=diffm, se, steps, bk)
  if (n<nmin) n <- nmin
  df  <- eval(dfe)
  pow <- .power.noninf(alpha=alpha, lmargin=lmargin, diffm=diffm, 
                       sem=se*sqrt(bk/n), df=df)
  if (details){
    cat("\nSample size search (ntotal)\n")
    cat(" n     power\n")
    # do not print first too high
    # this is for cases with only one step-down and than step up
    if (pow<=targetpower) cat( n," ", formatC(pow, digits=6, format="f"),"\n")
  }
  iter <- 0
  while (pow>targetpower){
    if (n<=nmin) {
      if (details & iter==0) cat( n," ", formatC(pow, digits=6, format="f"),"\n")
      break
    }
    n <- n - steps
    df  <- eval(dfe)
    pow <- .power.noninf(alpha=alpha, lmargin=lmargin, diffm=diffm, 
                         sem=se*sqrt(bk/n), df=df)
    iter <- iter+1
    if (details) cat( n," ", formatC(pow, digits=6, format="f"),"\n")
    if (iter>imax) break  
  }
  while(pow<targetpower){
    n <- n+steps
    df  <- eval(dfe)
    pow <- .power.noninf(alpha=alpha, lmargin=lmargin, diffm=diffm, 
                         sem=se*sqrt(bk/n), df=df)
    if (details) cat( n," ", formatC(pow, digits=6, format="f"),"\n")
    iter <- iter +1
    if (iter>imax) break  
  }
  if (pow<targetpower) {
    n <- NA
    if (details) cat("Sample size search failed!\n")
  }
  if (print && !details) {
    cat("\nSample size (total)\n")
    cat(" n     power\n")
    cat( n," ", formatC(pow, digits=6, format="f"),"\n")
    if (is.na(n)) cat("Sample size search failed!\n")
  }
  # return value: a data.frame
  res <- data.frame(design=design, alpha=alpha, CV=CV, theta0=theta0, 
                    margin=margin, n=n, power=pow, targetpower=targetpower)
  names(res) <- c("Design","alpha","CV","theta0","Margin", "Sample size", 
                  "Achieved power", "Target power")
  
  if (print) return(invisible(res)) else return(res)
}


