#---------------------------------------------------------------------------
# Sample size for partial and full replicate design with BE decision
# via ABE + sratio <=2.5 according to the FDA Dabigatran / rivaroxaban guidances
# via simulated (empirical) power
# 
# Author: dlabes
#---------------------------------------------------------------------------

sampleN.HVNTID <- function(alpha=0.05, targetpower=0.8, theta0, theta1, theta2, 
                           CV, design=c("2x2x4", "2x2x3"), nsims=1E5, nstart, 
                           imax=100, print=TRUE, details=TRUE, setseed=TRUE)
{
  if (missing(theta1) & missing(theta2)) theta1 <- 0.8
  if (missing(theta0)) theta0 <- 0.95
  if (missing(theta2)) theta2=1/theta1
  if ( (theta0<=theta1) | (theta0>=theta2) ) {
    stop("Null ratio ",theta0," not between margins ",theta1," / ",theta2,"!", 
         call.=FALSE)
  }
  if (missing(CV)) stop("CV(s) must be given!", call.=FALSE)
  
  # CV's for Test/Reference
  CVwT <- CV[1]
  if (length(CV)>1) CVwR <- CV[2] else CVwR <- CVwT
  if (length(CV)>2) warning("Only first 2 entries from CV vector used.")
  s2wT <- CV2mse(CVwT)
  s2wR <- CV2mse(CVwR)
  # check design argument
  design <- match.arg(design)
  # we are treating only balanced designs
  # thus we use here bk - the design constant for ntotal
  # expressions for the df's
  if (design=="2x2x4") {
    seqs  <- 2
    bk    <- 1.0    # needed for n0
    dfe   <- parse(text="n-2", srcfile=NULL)
    dfRRe <- parse(text="n-2", srcfile=NULL)
    dfTTe <- parse(text="n-2", srcfile=NULL)
    # expectation of mse of the ANOVA of intra-subject contrasts T-R
    Emse  <- (s2wT + s2wR)/2 
  }
  if (design=="2x2x3") {
    seqs  <- 2
    bk    <- 1.5    # needed for n0
    dfe   <- parse(text="n-2", srcfile=NULL)
    dfRRe <- parse(text="n/2-1", srcfile=NULL)  # balanced only
    dfTTe <- parse(text="n/2-1", srcfile=NULL)  # balanced only
    # expectation of mse of the ANOVA of intra-subject contrasts T-R
    Emse  <- 1.5*(s2wT + s2wR)/2                # balanced only
  }
  
  mlog <- log(theta0)
  ltheta1 <- log(theta1)
  ltheta2 <- log(theta2)
  # we have eventually to adapt this empirically for the n0 start value
  bkk     <- bk
  
  if (print){
    cat("\n+++++++++ FDA method for HV NTID's +++++++++++\n")
    cat("           Sample size estimation\n")
    cat("----------------------------------------------\n")
    cat("Study design: ",design,"\n")
    cat("log-transformed data (multiplicative model)\n")
    cat(nsims,"studies for each step simulated.\n\n")
    cat("alpha  = ",alpha,", target power = ", targetpower,"\n", sep="")
    cat("CVw(T) = ",CVwT,", CVw(R) = ",CVwR,"\n", sep="")
    cat("Null (true) ratio =",theta0,"\n")
    cat("ABE limits        =", theta1, "...", theta2,"\n")
#     if (details) {
#       cat("Implied scABEL    =", formatC(exp(ltheta1), format="f", digits=4), 
#           "...", formatC(exp(ltheta2), format="f", digits=4),"\n")
#     }
#     cat("Regulatory settings:",regulator,"\n")
#     if (details) { 
#       cat("- Regulatory const. =",r_const,"\n")
#       # CVcap?
#       cat("- 'CVcap'           =", formatC(CVcap, format="f", digits=4),"\n")
#     }     
  }
  # -----------------------------------------------------------------
  # nstart? from sampleN0 with shrunken/widened limits
  # does'nt fit always really good 
  if (missing(nstart)){
    n <- .sampleN0(alpha=alpha, targetpower, ltheta1, ltheta2, diffm=mlog, 
                   se=sqrt(Emse), steps=seqs, bk=bkk, diffmthreshold=0.01)
    #cat("n0=",n,"\n")
  } else n <- seqs*round(nstart/seqs)
  nmin <- 6 
  if(n<nmin) n <- nmin
  # we are simulating for balanced designs
  C3 <- 1/n
  # sd of the sample mean T-R (point estimator)
  sdm  <- sqrt(Emse*C3)
  df   <- eval(dfe)
  dfRR <- eval(dfRRe)
  dfTT <- eval(dfTTe)
  
  if(setseed) set.seed(123456)
  p <- .power.HVNTID(mlog, sdm, C3, Emse, df, s2wR, dfRR, s2wT, dfTT, nsims, 
                     ln_lBEL=log(theta1),ln_uBEL=log(theta2), alpha)
  pwr <- as.numeric(p["BE"]);
  
  if (details) {
    cat("\nSample size search\n")
    cat(" n     power\n")
    # do not print first too high
    # this is for cases with only one step-down and than step up
    if (pwr<=targetpower) cat( n," ", formatC(pwr, digits=6, format="f"),"\n")
  }
  iter <- 0
  # iter>100 is emergency brake
  # --- loop until power <= target power, step-down
  down <- FALSE; up <- FALSE
  while (pwr>targetpower) {
    if (n<=nmin) { 
      if (details & iter==0) cat( n," ", formatC(pwr, digits=6, format="f"),"\n")
      break
    }
    down <- TRUE
    n  <- n-seqs     # step down if start power is to high
    iter <- iter + 1
    C3 <- 1/n
    # sd of the sample mean T-R (point estimator)
    sdm  <- sqrt(Emse*C3)
    df   <- eval(dfe)
    dfRR <- eval(dfRRe)
    dfTT <- eval(dfTTe)
    
    if(setseed) set.seed(123456)
    p <- .power.HVNTID(mlog, sdm, C3, Emse, df, s2wR, dfRR, s2wT, dfTT, nsims, 
                       ln_lBEL=log(theta1),ln_uBEL=log(theta2), alpha)
    pwr <- as.numeric(p["BE"]);
    
    # do not print first step down
    if (details) cat( n," ", formatC(pwr, digits=6, format="f"),"\n")
    if (iter>imax) break  
    # loop results in n with power too low
    # must step one up again. is done in the next loop
  }
  while (pwr<targetpower) {
    up   <- TRUE ; down <- FALSE
    n    <- n+seqs  # step up
    iter <- iter+1
    C3 <- 1/n
    sdm  <- sqrt(Emse*C3)
    df   <- eval(dfe)
    dfRR <- eval(dfRRe)
    dfTT <- eval(dfTTe)
    
    if(setseed) set.seed(123456)
    p <- .power.HVNTID(mlog, sdm, C3, Emse, df, s2wR, dfRR, s2wT, dfTT, nsims, 
                       ln_lBEL=log(theta1),ln_uBEL=log(theta2), alpha)
    pwr <- as.numeric(p["BE"]);
    
    if (details) cat( n," ", formatC(pwr, digits=6, format="f"),"\n")
    if (iter>imax) break 
  }
  nlast <- n
  if (up & pwr<targetpower) {
    n <- NA
    if (details) cat("Sample size search failed!\n")
  }
  if (down & pwr>targetpower & n != nmin) {
    n <- NA
    if (details) cat("Sample size search failed!\n")
  }
  
  if (print && !details) {
    cat("\nSample size\n")
    cat(" n     power\n")
    cat( n," ", formatC(pwr, digits=6, format="f"),"\n")
    if (is.na(n)) cat("Sample size search failed!\n")
  }
  if (print) cat("\n")
  
  #return results as data.frame
  res <- data.frame(design=design, alpha=alpha, CVwT=CVwT, CVwR=CVwR,
                    theta0=theta0, theta1=theta1, theta2=theta2, n=n, power=pwr, 
                    targetpower=targetpower,nlast=nlast)
  names(res) <-c("Design","alpha","CVwT","CVwR","theta0","theta1","theta2",
                 "Sample size", "Achieved power", "Target power","nlast")
  
  if (print | details) return(invisible(res)) else return(res)
  
} # end function


