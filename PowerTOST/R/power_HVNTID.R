#---------------------------------------------------------------------------
# Simulate full replicate designs and calculate power of
# ABE + variability ratio <=2.5
# according to the FDA Dabigatran / rivaroxaban guidances
#
# Author: dlabes
#---------------------------------------------------------------------------

power.HVNTID <- function(alpha=0.05, theta1, theta2, theta0, CV, n, 
                         design=c("2x2x4", "2x2x3"), nsims=1E5, details=FALSE, 
                         setseed=TRUE)
{
  if (missing(CV)) stop("CV must be given!", call.=FALSE)
  if (missing(n))  stop("Number of subjects n must be given!", call.=FALSE)
   
  if (missing(theta0)) theta0 <- 0.95  
  if (missing(theta1) & missing(theta2)) theta1 <- 0.8
  if (missing(theta2)) theta2 <- 1/theta1
  # check design arg
  design <- match.arg(design)
  
  if(design=="2x2x4"){
    seqs  <- 2
    dfe   <- parse(text="n-2", srcfile=NULL)
    dfRRe <- parse(text="n-2", srcfile=NULL)
    dfTTe <- parse(text="n-2", srcfile=NULL)
  }
  if(design=="2x2x3"){
    seqs  <- 2
    dfe   <- parse(text="n-2", srcfile=NULL) 
    dfRRe <- parse(text="n/2-1", srcfile=NULL) # balanced only, not used here
    dfTTe <- parse(text="n/2-1", srcfile=NULL) # balanced only, not used here
  }
  
  CVwT <- CV[1]
  if (length(CV)>1) CVwR <- CV[2] else CVwR <- CVwT
  if (length(CV)>2) warning("Only first 2 entries from CV vector used.")
  s2wT <- CV2mse(CVwT)
  s2wR <- CV2mse(CVwR)

  if (length(n)==1){
    # we assume n=ntotal
    # for unbalanced designs we divide the ns by ourself
    # in such a way that we have only small imbalance
    nv <- nvec(n=n, grps=seqs)
    if (nv[1]!=nv[length(nv)]){
      message("Unbalanced design. n(i)=", paste(nv, collapse="/"), " assumed.")
    }
    C3 <- sum(1/nv)/seqs^2
    n  <- sum(nv)
  } else {
    # we assume n = vector of n's in sequence groups
    # check length
    if (length(n)!=seqs) stop("n must be a vector of length=",seqs,"!", call.=FALSE)
    nv <- n
    C3 <- sum(1/n)/seqs^2
    n  <- sum(n)
  }
  
  if(design=="2x2x4"){
    dfRR <- eval(dfRRe)
    dfTT <- dfRR       
    # expectation of mse of the ANOVA of intra-subject contrasts T-R
    Emse  <- (s2wT + s2wR)/2 
  }
  if(design=="2x2x3"){
    dfTT <- nv[1]-1
    dfRR <- nv[2]-1
    w1 <- dfRR/(dfRR+dfTT); w2 <- dfTT/(dfRR+dfTT)
    # expectation of mse of the ANOVA of intra-subject contrasts T-R
    # always via unbalanced formula
    Emse <- w1*(s2wT+s2wR/2) + w2*(s2wT/2+s2wR)
  }
  
  # start time measurement
  ptm <- proc.time()
  
  df <- eval(dfe)
  # sd of the mean T-R (point estimator)
  sdm  <- sqrt(Emse*C3)
  mlog <- log(theta0)
  
  if(setseed) set.seed(123456)
  
  p <- .power.HVNTID(mlog, sdm, C3, Emse, df, s2wR, dfRR, s2wT, dfTT, nsims, 
                     ln_lBEL=log(theta1),ln_uBEL=log(theta2), alpha=alpha)
    
  if (details) {
    ptm <- summary(proc.time()-ptm)
    message(nsims,"sims. Time elapsed (sec): ", 
            formatC(ptm["elapsed"], digits=2), "\n")
    #print(ptm)
    
    # return power and components
    names(p) <- c("p(BE)", "p(BE-ABE)", "p(BE-sratio)")
    p
  } else {
    # return only the 'power', overall p(BE)
    as.numeric(p["BE"])
  }
}

# working horse of RSABE for NTID's
.power.HVNTID <- function(mlog, sdm, C3, Emse, df, s2wR, dfRR, s2wT, dfTT, nsims, 
                          ln_lBEL=log(0.8), ln_uBEL=log(1.25), alpha=0.05)
{
  tval     <- qt(1-alpha,df)
  Fval     <- qf(1-alpha, dfTT, dfRR, lower.tail=FALSE)
  
  counts   <- rep.int(0, times=3)
  names(counts) <- c("BE", "BEabe", "BEsratio")
  # to avoid memory problems for high number of sims
  # we are working with chunks of 1e7
  chunks   <- 1
  nsi      <- nsims
  if (nsims>1E7) {
    chunks <- round(nsims/1E7,0)
    nsi    <- 1E7
  } 
  for (iter in 1:chunks) {
    # simulate sample mean via its normal distribution
    means  <- rnorm(nsi, mean=mlog, sd=sdm)
    # simulate sample sd2s via chi-square distri
    sd2s   <- Emse*C3*rchisq(nsi, df)/df
    # simulate sample value s2wRs via chi-square distri
    s2wRs  <- s2wR*rchisq(nsi, dfRR)/dfRR
    # simulate sample value s2wTs via chi-square distri
    s2wTs  <- s2wT*rchisq(nsi, dfTT)/dfTT
    
    SEs <- sqrt(sd2s)
    
    # conventional (1-2*alpha) CI's for T-R
    hw  <- tval*SEs
    lCL <- means - hw 
    uCL <- means + hw
    # conventional ABE
    BEABE     <- ((ln_lBEL<=lCL) & (uCL<=ln_uBEL))
    
    # save memory
    rm(SEs, hw)
    
    # upper limit of ratio swT/swR
    ul_sratio <- sqrt(s2wTs/s2wRs/Fval)
    # upper limit <= 2.5?
    BEsratio  <- ul_sratio <= 2.5
    
    counts["BEabe"]    <- counts["BEabe"]    + sum(BEABE)
    counts["BEsratio"] <- counts["BEsratio"] + sum(BEsratio)
    counts["BE"]       <- counts["BE"]       + sum(BEABE & BEsratio)
    
  } # end over chunks
  
  # return the pBEs
  counts/nsims
}