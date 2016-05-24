#---------------------------------------------------------------------------
# Sample size for partial and full replicate design and scaled ABE 
# via simulated (empirical) power
# 
# Author: dlabes
#---------------------------------------------------------------------------

sampleN.RSABE <- function(alpha=0.05, targetpower=0.8, theta0, theta1, 
                          theta2, CV, design=c("2x3x3", "2x2x4", "2x2x3"),
                          regulator = c("FDA", "EMA"), nsims=1E5, nstart, imax=100,
                          print=TRUE, details=TRUE, setseed=TRUE)
{
  if (missing(theta1) & missing(theta2)) theta1 <- 0.8
  # according to the two Laszlo's paper 0.9 (alas 1.10) should be considered
  if (missing(theta0)) theta0 <- 0.90
  if (missing(theta2)) theta2=1/theta1
  if ( (theta0<=theta1) | (theta0>=theta2) ) {
    stop("Null ratio ",theta0," not between margins ",theta1," / ",theta2,"!", 
         call.=FALSE)
  }
  if (missing(CV)) stop("CV(s) must be given!", call.=FALSE)
  
  #if (!print) details <- FALSE # do not print anything?

  CVswitch  <- 0.3
  regulator <- match.arg(regulator)
  if (regulator=="FDA") r_const <- log(1.25)/0.25 # or better log(theta2)/0.25?
  if (regulator=="EMA") r_const <- 0.76 # or better log(theta2)/CV2se(0.3)

  # for later enhancement taking into account the 
  # subject-by-formulation interaction
  s2D  <- 0 
  CVwT <- CV[1]
  if (length(CV)==2) CVwR <- CV[2] else CVwR <- CVwT
  s2wT <- log(1.0 + CVwT^2)
  s2wR <- log(1.0 + CVwR^2)
  
  # check design
  design <- match.arg(design)
  # we are treating only balanced designs
  # thus we use here bk - design constant for ntotal
  # expressions for the df's
  if (design=="2x3x3") {
    desi <- "2x3x3 (partial replicate)"
    seqs <- 3
    bk   <- 1.5    # needed for n0
    # in case of the FDA we are using the 'robust' df's
    # due to the fact that the described analysis in the
    # progesterone guidance is based in the intrasubject contrasts
    # T-R and R-R with df=n-seqs
    dfe   <- parse(text="n-3", srcfile=NULL)
    dfRRe <- parse(text="n-3", srcfile=NULL)
    # expectation of mse of the ANOVA of intra-subject contrasts
    #sd2  <- s2D + (s2wT + s2wR)/2 # used in v1.1-00 - v1.1-02
    # according to McNally et al., verified via simulations:
    Emse  <- s2D + s2wT + s2wR/2
  }
  if (design=="2x2x4") {
    desi <- "2x2x4 (full replicate)"
    seqs <- 2
    bk   <- 1.0    # needed for n0
    dfe   <- parse(text="n-2", srcfile=NULL)
    dfRRe <- parse(text="n-2", srcfile=NULL)
    # expectation of mse of the ANOVA of intra-subject contrasts
    Emse  <- (s2D + (s2wT + s2wR)/2) 
  }
  if (design=="2x2x3") {
    desi <- "2x2x3 (TRT|RTR)"
    seqs <- 2
    bk   <- 1.5    # needed for n0?
    dfe   <- parse(text="n-2", srcfile=NULL)
    dfRRe <- parse(text="n/2-1", srcfile=NULL) # for balanced designs
    # expectation of mse of the ANOVA of intra-subject contrasts
    Emse  <- 1.5*(s2wT + s2wR)/2               # for balanced designs
  }
  
  mlog <- log(theta0)
  
  if (print){
    cat("\n++++++++ Reference scaled ABE crit. +++++++++\n")
    cat("           Sample size estimation\n")
    cat("---------------------------------------------\n")
    cat("Study design: ",desi,"\n")
    cat("log-transformed data (multiplicative model)\n")
    cat(nsims,"studies for each step simulated.\n\n")
    cat("alpha  = ",alpha,", target power = ", targetpower,"\n", sep="")
    cat("CVw(T) = ",CVwT,"; CVw(R) = ",CVwR,"\n", sep="")
    cat("Null (true) ratio = ",theta0,"\n", sep="")
    cat("ABE limits / PE constraints =",theta1,"...", theta2,"\n")
    cat("Regulatory settings:",regulator,"\n")
    if (details) { 
      cat("- CVswitch = ", CVswitch, "\n")
      cat("- Regulatory constant =",r_const,"\n")
    }     
  }
  
  # -----------------------------------------------------------------
  # nstart from sampleN0 with widened limits
  # does'nt fit really good if theta0>=1.2 or <=0.85! ways out?
  ltheta1 <- -r_const*sqrt(s2wR)
  ltheta2 <- -ltheta1
  # this does not function in case of CVwR=0.3 for the original code
  # calculating s2wR and back-calculating CVwR from that
  # numerical problem?
  if (CVwR <= CVswitch){
    ltheta1 <- log(theta1)
    ltheta2 <- log(theta2)
  }
  if (missing(nstart)){
    
    # try to use the empirical alpha for start sample size, some sort of
    al <- alpha
    if (regulator=="FDA") {
      if(Emse/bk <= CV2mse(0.30001) & Emse/bk >= CV2mse(0.2975)) al=0.12
      if(Emse/bk > CV2mse(0.30001)) al <- 0.035   
    }
    if (regulator=="EMA") {
      #does not fit!
      #if(Emse/bk <= CV2mse(0.321) & Emse/bk >= CV2mse(0.28)) al=0.065
    }
    # debug print
    # cat(al,"\n")
    # we use bk=1 here since our formula is sem=sqrt(Emse/n)
    n01 <- .sampleN0(alpha=al, targetpower, ltheta1, ltheta2, diffm=mlog, 
                     se=sqrt(Emse), steps=seqs, bk=1, diffmthreshold=0.01)
    # empirical correction in the vicinity of CV=0.3, for ratios 
    # outside 0.86 ... 1/0.86
#     if(Emse/bk <= CV2mse(0.305) & Emse/bk >= CV2mse(0.295) & abs(mlog)>log(1/0.865)) {
#       if (regulator=="EMA") n01 <- 0.9*n01 else  n01 <- 0.8*n01
#       n01 <- seqs*trunc(n01/seqs)
#     }  
    # start from PE constraint sample size
    n02 <- .sampleN0.2(targetpower, ltheta2=log(theta2), diffm=mlog, 
                       se=sqrt(Emse), steps=seqs, bk=1)
    # debug print
    # cat(n01,n02,"\n")
    n <- max(c(n01,n02))+seqs
  } else n <- seqs*round(nstart/seqs)
  nmin <- 6 # fits 2x3x3 and 2x2x4
  if(n<nmin) n <- nmin
  # we are simulating for balanced designs
  C3 <- 1/n
  # sd of the sample mean T-R (point estimator)
  sdm  <- sqrt(Emse*C3)
  df   <- eval(dfe)
  dfRR <- eval(dfRRe)
  
  if(setseed) set.seed(123456)
  p <- .power.RSABE(mlog, sdm, C3, Emse, df, s2wR, dfRR, nsims, 
                    ln_lBEL=log(theta1), ln_uBEL=log(theta2), 
                    CVswitch, r_const, alpha=alpha)
  pwr <- as.numeric(p["BE"]);
  pd <- max(4,round(log10(nsims),0))  # digits for power
  if (details) {
    cat("\nSample size search\n")
    cat(" n     power\n")
    # do not print first too high
    # this is for cases with only one step-down and than step up
    if (pwr<=targetpower) cat( n," ", formatC(pwr, digits=pd, format="f"),"\n")
  }
  iter <- 0
  # iter>100 is emergency brake
  # --- loop until power <= target power, step-down
  down <- FALSE; up <- FALSE
  while (pwr>targetpower) {
    down <- TRUE
    if (n<=nmin) { 
      if (details & iter==0) cat(n," ", formatC(pwr, digits=pd, format="f"),"\n")
      break
    }
    n  <- n-seqs     # step down if start power is to high
    iter <- iter + 1
    C3 <- 1/n
    # sd of the sample mean T-R (point estimator)
    sdm  <- sqrt(Emse*C3)
    df   <- eval(dfe)
    dfRR <- eval(dfRRe)
    
    if(setseed) set.seed(123456)
    p <- .power.RSABE(mlog, sdm, C3, Emse, df, s2wR, dfRR, nsims, 
                      ln_lBEL=log(theta1),ln_uBEL=log(theta2), 
                      CVswitch, r_const, alpha=alpha)
    pwr <- as.numeric(p["BE"]);
    
    # do not print first step down
    if (details) cat( n," ", formatC(pwr, digits=pd, format="f"),"\n")
    if (iter>imax) break  
    # loop results in n with power too low
    # must step one up again. is done in the next loop
  }
  while (pwr<targetpower) {
    down <- FALSE; up <- TRUE
    n    <- n+seqs  # step up
    iter <- iter+1
    C3 <- 1/n
    sdm  <- sqrt(Emse*C3)
    df   <- eval(dfe)
    dfRR <- eval(dfRRe)
    
    if(setseed) set.seed(123456)
    p <- .power.RSABE(mlog, sdm, C3, Emse, df, s2wR, dfRR, nsims, 
                      ln_lBEL=log(theta1),ln_uBEL=log(theta2), 
                      CVswitch, r_const, alpha=alpha)
    pwr <- as.numeric(p["BE"]);
    
    if (details) cat( n," ", formatC(pwr, digits=pd, format="f"),"\n")
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
    cat(" n    power\n")
    cat( n," ", formatC(pwr, digits=pd, format="f"),"\n")
    if (is.na(n)) cat("Sample size search failed!\n")
  }
  if (print) cat("\n")
  
  #return results as data.frame
  res <- data.frame(design=design, alpha=alpha, CVwT=CVwT, CVwR=CVwR,
                    theta0=theta0, theta1=theta1, theta2=theta2, n=n, power=pwr, 
                    targetpower=targetpower,nlast=nlast)
  names(res) <-c("Design","alpha","CVwT","CVwR","theta0","theta1","theta2",
                 "Sample size", "Achieved power", "Target power","nlast")
  # debug print
  # cat("iter=",iter+1,"\n")
  
  if (print | details) return(invisible(res)) else return(res)
} # end function


