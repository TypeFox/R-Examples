#-----------------------------------------------------------------------------
# Author: dlabes
#-----------------------------------------------------------------------------

# ----- helper functions for sampleN.TOST and other --------------------------
# Sample size for a desired power, large sample approx.
# bk = design constant, see known.designs()
.sampleN0 <- function(alpha=0.05, targetpower=0.8, ltheta1, ltheta2, diffm, 
                      se, steps=2, bk=2, diffmthreshold=0.03)
{
  z1 <- qnorm(1-alpha)
  # value diffmthreshold=0.04 corresponds roughly to log(0.96)
  # with lower values there are many steps around between 0.95 and 1
  # in sampleN.TOST
  if (abs(diffm)>diffmthreshold) z2 <- qnorm(targetpower) else {
    z2 <- qnorm(1-(1-targetpower)*0.5) # for diffm ~0 (log: theta0=1) 1-beta/2
    diffm <- 0
  }
  n01<-(bk/2)*((z1+z2)*(se*sqrt(2)/(diffm-ltheta1)))^2;
  n02<-(bk/2)*((z1+z2)*(se*sqrt(2)/(diffm-ltheta2)))^2;
  
  n0 <- ceiling(pmax(n01,n02))
  #make an even multiple of step (=2 in case of 2x2 cross-over)
  n0 <- steps*trunc(n0/steps)
  
  # minimum sample size will be checked outside
  return(n0)
}

# pure Zhang's formula. doesn't work sufficiently?
# Why? Don't know any more DL Jan 2015
# examine_n0 shows superiority of this function over next one for 2x2
# Paul Zhang (2003)
# A Simple Formula for Sample Size Calculation in Equivalence Studies 
# Journal of Biopharmaceutical Statistics, 13:3, 529-538
.sampleN0_2 <- function(alpha=0.05, targetpower=0.8, ltheta1, ltheta2, diffm, 
                        se, steps=2, bk=2)
{
  # handle unsymmetric limits, Zhang's c0
  c0 <- 0.5*exp(-7.06*(ltheta1+ltheta2)/(ltheta1-ltheta2))
  # Zhang's formula, large sample
  beta <- 1-targetpower
  z1 <- qnorm(1-alpha)
  fz <- ifelse(diffm<0, c0*exp(-7.06*diffm/ltheta1), c0*exp(-7.06*diffm/ltheta2))
  z2 <- abs(qnorm((1-fz)*beta))
  
  n01<-(bk/2)*((z1+z2)*(se*sqrt(2)/(diffm-ltheta1)))^2;
  n02<-(bk/2)*((z1+z2)*(se*sqrt(2)/(diffm-ltheta2)))^2;
  
  n0 <- pmax(n01,n02) # or ceiling/round?
  # make an even multiple of step (=2 in case of 2x2 cross-over)
  n0 <- steps*trunc(n0/steps)
  
  return(n0)
}

# mixture of old code and Zhang's formula 
.sampleN0_3 <- function(alpha=0.05, targetpower=0.8, ltheta1, ltheta2, diffm, 
                        se, steps=2, bk=2)
{
  # transform to limits symmetric around zero (if they are not)
  locc    <- (ltheta1+ltheta2)/2
  diffm   <- diffm - locc
  ltheta1 <- ltheta1 - locc
  ltheta2 <- -ltheta1
  delta   <- abs((ltheta2-ltheta1)/2)
  
  z1   <- qnorm(1-alpha)
  beta <- 1-targetpower
  
  c  <- abs(diffm/delta)
  # probability for second normal quantil
  # c=0.2 corresponds roughly to exp(diffm) >0.95 or < 1.05 if ltheta1/ltheta2
  # are the logs of 0.8/1.25
  # c=0.35 corresponds roughly to >0.925 ... < 1.08
  # outside these we use 1-beta, inside smooth change to 1-beta/2
  p2 <- ifelse(c<0.35, 1-(1-0.5*exp(-7.06*c))*beta, 1-beta)
  z2 <- qnorm(p2)
  # difference for denominator
  dn <- ifelse(diffm<0, diffm-ltheta1, diffm-ltheta2)
  n0 <- (bk/2)*((z1+z2)*(se*sqrt(2)/dn))^2
  # make an even multiple of steps (=2 in case of 2x2 cross-over)
  n0 <- steps*trunc(n0/steps)
  
  return(n0)
}

#------------------------------------------------------------------------------
# Sample size for a desired power: 
# see known.designs() for covered experimental designs
# theta1 if empty is set to 0.8 or -0.2 depending on logscale
# diff if empty is set to 0.95 or 0.05 depending on logscale
# leave upper BE margin (theta2) empty and the function will use -lower
# in case of additive model or 1/lower if logscale=TRUE
sampleN.TOST <- function(alpha=0.05, targetpower=0.8, logscale=TRUE, theta0, 
                         theta1, theta2, CV, design="2x2", method="exact",
                         robust=FALSE, print=TRUE, details=FALSE, imax=100)
{
  if (missing(CV)) stop("CV must be given!", call.=FALSE)
  if(any(CV<0)) {
    message("Negative CV changed to abs(CV).")
    CV <- abs(CV)
  }
  
  #number of the design and check
  d.no <- .design.no(design)
  if (is.na(d.no)) stop("Design ",design, " unknown!", call.=FALSE)
  
  # design characteristics
  ades   <- .design.props(d.no)
  d.name <- ades$name  # nice name of design
  # get the df for the design as an unevaluated expression (now with n as var)
  dfe    <- .design.df(ades, robust=robust)
  steps  <- ades$steps	# stepsize for sample size search
  bk     <- ades$bk     # get design constant
  # minimum n
  df <- 0
  n  <- 0
  while (df<1){
    n  <- n + 1
    df <- eval(dfe)
  }
  # make a multiple of steps
  nmin <- as.integer(steps*trunc(n/steps)) 
  nmin <- nmin + steps*(nmin<n)
  # print the configuration:
  if (print) {
    cat("\n+++++++++++ Equivalence test - TOST +++++++++++\n")
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
  
  # regularize the method giving
  method <- .powerMethod(method)
  
  # handle the log transformation
  if (logscale) {
    if (missing(theta1) & missing(theta2)) theta1 <- 0.8
    if (missing(theta0)) theta0 <- 0.95
    if (missing(theta2)) theta2=1/theta1
    if ( (theta0<=theta1) | (theta0>=theta2) ) {
      stop("Null ratio ",theta0," not between margins ",theta1," / ",theta2,"!", 
           call.=FALSE)
     }
    ltheta1 <- log(theta1)
    ltheta2 <- log(theta2)
    diffm   <- log(theta0)
    se      <- CV2se(CV)
    if (print) cat("log-transformed data (multiplicative model)\n\n")
  } else {
    if (missing(theta1) & missing(theta2)) theta1 <- -0.2
    if (missing(theta0)) theta0 <- 0.05
    if (missing(theta2)) theta2=-theta1
    if ( (theta0<=theta1) | (theta0>=theta2) ) {
      stop("Null diff. ",theta0," not between margins ",theta1," / ",theta2,"!", 
           call.=FALSE)
    }
    ltheta1 <- theta1
    ltheta2 <- theta2
    diffm   <- theta0
    se      <- CV
    if (print) cat("untransformed data (additive model)\n\n")
  }
  
  if (print) {
    cat("alpha = ",alpha,", target power = ", targetpower,"\n", sep="")
    cat("BE margins        =",theta1,"...", theta2,"\n")
    if (logscale) cat("Null (true) ratio = ",theta0,",  CV = ",CV,"\n", sep="")
    else          cat("Null (true) diff. = ",theta0,",  CV = ",CV,"\n", sep="")
  }
  
  # start value from large sample approx. (hidden func.)
  # Jan 2015 attempt to change to pure Zhang's formula
  # but this gives many iterations if theta0 near acceptance limits
  n <- .sampleN0_3(alpha, targetpower, ltheta1, ltheta2, diffm, se, steps, bk)
  if (n<nmin) n <- nmin

  df <- eval(dfe)
  pow <- .calc.power(alpha=alpha, ltheta1=ltheta1, ltheta2=ltheta2, diffm=diffm, 
                     sem=se*sqrt(bk/n), df=df, method=method)

  if (details) {
    cat("\nSample size search (ntotal)\n")
    cat(" n     power\n")
    # do not print first too high
    # this is for cases with only one step-down and than step up
    if (pow<=targetpower) cat( n," ", formatC(pow, digits=6, format="f"),"\n")
  }
  iter <- 0
  # iter>imax is emergency brake
  # this is eventually not necessary, depends on quality of sampleN0
  # in experimentation I have seen max of six steps
  # reformulation with only one loop does not shorten the code considerable
  # --- loop until power <= target power, step-down
  down <- FALSE; up <- FALSE
  while (pow>targetpower) {
    if (n<=nmin) { 
      if (details & iter==0) cat( n," ", formatC(pow, digits=6, format="f"),"\n")
      break
    }
    down <- TRUE
    n    <- n-steps     # step down if start power is to high
    iter <- iter+1
    df   <- eval(dfe)
    pow <- .calc.power(alpha=alpha, ltheta1=ltheta1, ltheta2=ltheta2, diffm=diffm, 
                       sem=se*sqrt(bk/n), df=df, method=method)
      
    # do not print first step down
    if (details) cat( n," ", formatC(pow, digits=6),"\n")
    if (iter>imax) break  
    # loop results in n with power too low
    # must step one up again. is done in the next loop
  }
  # --- loop until power >= target power
  while (pow<targetpower) {
    up   <- TRUE; down <- FALSE
    n    <- n+steps
    iter <- iter+1
    df   <- eval(dfe)
    pow <- .calc.power(alpha, ltheta1, ltheta2, diffm, sem=se*sqrt(bk/n), df, method)
    if (details) cat( n," ", formatC(pow, digits=6, format="f"),"\n")
    if (iter>imax) break 
  }
  #browser()
  # for very large n the search usually fails, but doesn't matter
  nlast <- n
  if (up & pow<targetpower & n<10000) {
    n <- NA
    if (details) cat("Sample size search failed!\n")
  }
  if (down & pow>targetpower & n<10000) {
    n <- NA
    if (details) cat("Sample size search failed!\n")
  } 
  if (print && !details) {
    cat("\nSample size (total)\n")
    #if (d.no == 0) cat("(n is sample size per group)\n") #parallel group design
    cat(" n     power\n")
    cat( n," ", formatC(pow, digits=6, format="f"),"\n")
    if (is.na(n)) cat("Sample size search failed!\n")
  }
  if (details && print) {
    if (method=="exact" || method=="owenq") 
      cat("\nExact power calculation with\nOwen's Q functions.\n")
    if (method=="mvt") 
      cat("\nExact power calculation with\nbivariate non-central t-distribution.\n")
  }
  # always print if approx.
  if (print & (method!="exact" & method!="mvt")){
    approx <- switch(
      method,
      nct="Approximate power calculation with\nnon-central t-distribution.",
      noncentral="Approximate power calculation with\nnon-central t-distribution.",
      shifted="Approximate power calculation with\nshifted central t-distribution.",
      central="Approximate power calculation with\nshifted central t-distribution."
      )
    cat("\n",approx,"\n",sep="")
  } 
  if (print) cat("\n")
  
  #return results as data.frame
  res <- data.frame(design=design, alpha=alpha, CV=CV, theta0=theta0, 
                    theta1=theta1, theta2=theta2, n=n, power=pow, 
                    targetpower=targetpower)
  names(res) <-c("Design","alpha","CV","theta0","theta1","theta2",
                 "Sample size", "Achieved power", "Target power")
  
  if (print) return(invisible(res)) else return(res)
  
}
