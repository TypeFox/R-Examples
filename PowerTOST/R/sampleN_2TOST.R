#-----------------------------------------------------------------------------
# Author: dlabes
# Adapted for power.2TOST by Benjamin Lang
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
sampleN.2TOST <- function(alpha=c(0.05, 0.05), targetpower=0.8, logscale=TRUE, 
                          theta0, theta1, theta2, CV, rho, design="2x2", 
                          setseed=TRUE, robust=FALSE, print=TRUE, details=FALSE,
                          imax=100)
{
  if (length(alpha) != 2)
    stop("Two alpha values must be given!")
  if (!missing(theta0) && length(theta0) != 2)
    stop("Two theta0 values must be given!")
  if (!missing(theta1) && length(theta1) != 2)
    stop("Two theta1 values must be given!")
  if (!missing(theta2) && length(theta2) != 2)
    stop("Two theta2 values must be given!")
  if (missing(CV)) stop("CV must be given!", call.=FALSE)
  if (length(CV) != 2)
    stop("Two CVs must be given!")
  if(any(CV<0)) {
    message("Negative CV changed to abs(CV).")
    CV <- abs(CV)
  }
  if (missing(rho))
    stop("Correlation between the two endpoints must be given!")
  if (length(rho) != 1)
    stop("One rho must be given!")
  if (rho < -1 || rho > 1)
    stop("Correlation must be >= -1 and =< 1.") 
  
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
    cat("\n+++++++++++ Equivalence test - 2 TOSTs +++++++++++\n")
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
    if (missing(theta1)) theta1 <- c(0.8, 0.8)
    if (missing(theta2)) theta2 <- 1/theta1
    if (any(theta1 <= 0) || any(theta1 > theta2))
      stop("theta1 and/or theta2 not correctly specified.")
    if (missing(theta0)) theta0 <- c(0.95, 0.95)
    if (any(theta0 <= 0))
      stop("theta0 must be > 0.")
    if (any(theta0 <= theta1) || any(theta0 >= theta2)) {
      stop("One assumed ratio is not between respective margins!", 
           call. = FALSE)
    }
    ltheta1 <- log(theta1)
    ltheta2 <- log(theta2)
    diffm   <- log(theta0)
    se      <- CV2se(CV)
    if (print) cat("log-transformed data (multiplicative model)\n\n")
  } else {
    if (missing(theta0)) 
      theta0 <- c(0.05, 0.05)
    if (missing(theta1)) 
      theta1 <- c(-0.2, -0.2)
    if (missing(theta2)) 
      theta2 <- -theta1
    if (any(theta1 > theta2))
      stop("theta1 and/or theta2 not correctly specified.")
    if (any(theta0 <= theta1) || any(theta0 >= theta2)) {
      stop("One assumed difference is not between respective margins!", 
           call. = FALSE)
    }
    ltheta1 <- theta1
    ltheta2 <- theta2
    diffm   <- theta0
    se      <- CV
    if (print) cat("untransformed data (additive model)\n\n")
  }
  
  if (print) {
    cat("alpha = ",paste(as.character(alpha), collapse = ", "),
        "; target power = ", targetpower,"\n", sep="")
    cat("BE margins         =",theta1,"...", theta2,"\n")
    if (logscale) cat("Null (true) ratios = ",paste(as.character(theta0), 
                                                    collapse = ", "),
                      "; CV = ",paste(as.character(CV), collapse = ", "),
                      "\n", sep="")
    else          cat("Null (true) diffs = ",paste(as.character(theta0), 
                                                    collapse = ", "),
                      "; SD = ",paste(as.character(CV), collapse = ", "),
                      "\n", sep="")
    cat("Correlation between the two parameters = ",rho,"\n", sep="")
  }
  
  # if both theta0 are near acceptance limits then starting value may not be
  # ideal resulting in a lot of iteration steps
  idx.d <- which.max(abs(diffm))
  n <- n.tmp <- max(.sampleN0_3(alpha[1], targetpower, ltheta1[1], 
                                ltheta2[1], diffm[1], se[1], steps, bk),
                    .sampleN0_3(alpha[2], targetpower, ltheta1[2], 
                                ltheta2[2], diffm[2], se[2], steps, bk))
  df <- eval(dfe)
  pow <- .prob.2TOST(ltheta0 = diffm, se = se*sqrt(bk/n), df = df,
                      ltheta1 = ltheta1, ltheta2 = ltheta2, rho = rho,
                      alpha = alpha, setseed = setseed)
  if (!isTRUE(all.equal(pow, targetpower, tolerance = 1e-04))) {
    n <- .sampleN0_3(min(alpha), targetpower, ltheta1[idx.d], ltheta2[idx.d], 
                     diffm[idx.d], max(se), steps, bk)
    df <- eval(dfe)
    pow.tmp <- .prob.2TOST(ltheta0 = diffm, se = se*sqrt(bk/n), df = df,
                            ltheta1 = ltheta1, ltheta2 = ltheta2, rho = rho,
                            alpha = alpha, setseed = setseed)
    if (abs(pow.tmp - targetpower) <= abs(pow - targetpower)) {
      pow <- pow.tmp
    } else {
      n <- n.tmp
      df <- eval(dfe)
    }
  }
  if (n < nmin) n <- nmin
  
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
    n    <- n-steps     # step down if start power is too high
    iter <- iter+1
    df   <- eval(dfe)
    pow <- .prob.2TOST(ltheta0 = diffm, se = se*sqrt(bk/n), df = df,
                        ltheta1 = ltheta1, ltheta2 = ltheta2, rho = rho,
                        alpha = alpha, setseed = setseed)
    # do not print first step down
    if (details) cat( n," ", formatC(pow, digits=6),"\n")
    if (iter>imax) break  
    # loop results in n with power too low
    # must step one up again. is done in the next loop
  }
  # --- loop until power >= target power
  while (pow<targetpower) {
    up   <- TRUE; down <- FALSE
    n <- n+steps
    iter <- iter+1
    df   <- eval(dfe)
    pow <- .prob.2TOST(ltheta0 = diffm, se = se*sqrt(bk/n), df = df,
                        ltheta1 = ltheta1, ltheta2 = ltheta2, rho = rho,
                        alpha = alpha, setseed = setseed)
    if (details) cat( n," ", formatC(pow, digits=6, format="f"),"\n")
    if (iter>imax) break 
  }
  
  nlast <- n
  if (up & pow<targetpower) {
    n <- if (isTRUE(all.equal(pow, targetpower, 1e-04))) nlast else NA
    if (details && is.na(n)) 
      cat("Sample size search failed; last n = ", nlast,
          " (power = ", pow, ")\n", sep="")
  }
  if (down & pow>targetpower) {
    n <- if (isTRUE(all.equal(pow, targetpower, 1e-04))) nlast else NA
    if (details && is.na(n)) 
      cat("Sample size search failed; last n = ", nlast,
          " (power = ", pow, ")\n", sep="")
  } 
  if (print && !details) {
    cat("\nSample size (total)\n")
    #if (d.no == 0) cat("(n is sample size per group)\n") #parallel group design
    cat(" n     power\n")
    cat( n," ", formatC(pow, digits=6, format="f"),"\n")
    if (is.na(n)) 
      cat("Sample size search failed; last n = ", nlast, 
          " (power = ", pow, ")\n", sep="")
  }
  if (print) cat("\n")
  
  res <- list(design=design, alpha=alpha, CV=CV, theta0=theta0, theta1=theta1,
              theta2=theta2, n=n, power=pow, targetpower=targetpower)
  names(res) <- c("Design","alpha","CV","theta0","theta1","theta2",
                  "Sample size", "Achieved power", "Target power")
  
  if (print) return(invisible(res)) else return(res)
}