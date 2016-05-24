#---------------------------------------------------------------------------
# Simulate partial and full replicate design and pure scaled ABEL power
# 
# Author: dlabes
#---------------------------------------------------------------------------

# degrees of freedom for the RR  analysis: 
# model with subj, period like EMA Q&A
# the inclusion of sequence changes only the split of df's
# between seq with df=seq-1 and sub(seq) with df=(n-1) - (seq-1). 
# 2x3x3  2*n measurements df = 2*n-1
#        n subjects       df = n-1
#        3 periods        df = 2
#                   ->  dfRR = n-2
# 2x2x4  2*n measurements df = 2*n-1
#        n subjects       df = n-1
#        4 periods        df = 3     2
#                   ->  dfRR = n-3   n-2
# But cave! The EMA set I (2x2x4) has only 2 df for period. 
# Due to imbalance? No. Simulated data have also df=2 for period (typIII).
# for typI we have n-2 for subject and df=3 for period!
#
# Another possibility is using the contrasts R-R and analyze by sequence. 
# Then the dfRR = n-seq.
# 2x3x3  dfRR = n-3
# 2x2x2  dfRR = n/2 - 1 (only sequence TRR or RTR)
# 2x2x4  dfRR = n-2
# This is used in the xxx.RSABE() functions

pwr.scABEL <- function(alpha=0.05, theta1, theta2, theta0, CV, n,   
                       design=c("2x3x3", "2x2x4", "2x2x3"), 
                       regulator=c("EMA", "ANVISA", "FDA"), robust=FALSE,
                       nsims=1E5, setseed=TRUE)
{
  if (missing(CV)) stop("CV must be given!")
  if (missing(n))  stop("Number of subjects n must be given!")

  if (missing(theta0)) theta0 <- 0.95
  if (missing(theta1) & missing(theta2)) theta1 <- 0.8
  if (missing(theta2)) theta2 <- 1/theta1
  if (missing(theta1)) theta1 <- 1/theta2
  
  ptm <- proc.time()
  
  # subject-by-formulation interaction can't play a role here
  # since the model doesn't allow such term
  CVwT <- CV[1]
  if (length(CV)==2) CVwR <- CV[2] else CVwR <- CVwT
  s2wT <- log(1.0 + CVwT^2)
  s2wR <- log(1.0 + CVwR^2)
  
  regulator <- toupper(regulator)
  regulator <- match.arg(regulator)
  # constants acc. to regulatory bodies
  rc <- reg_const(regulator)
  CVcap    <- rc$CVcap
  CVswitch <- rc$CVswitch
  r_const  <- rc$r_const
  
  design <- match.arg(design)
  if (design=="2x3x3") {
    seqs <- 3
    bkni <- 1/6
    bk   <- 1.5
    dfe   <- parse(text="2*n-3", srcfile=NULL)
    if(robust) dfe <- parse(text="n-3", srcfile=NULL)
    dfRRe <- parse(text="n-2", srcfile=NULL)
    # sd^2 (variance) of the differences T-R from their components
    # According to McNally (sD=0):
    #sd2  <- s2wT + s2wR/2
    # but this gives too low power compared to the tables of the 2 Laszlos
    #sd2  <- (s2wT + s2wR)/2 # used in v1.1-00 - v1.1-02, wrong
    # simulations with s2D=0 show:
    Emse  <- (s2wT + 2.0*s2wR)/3
    cvec <- c(1,2)
    # warning in case of CVwR != CVwT
    if (abs(s2wT-s2wR)>1e-5){
      warning(paste("Heteroscedasticity in the 2x3x3 design may led", 
              "to poor accuracy of power!"), call.=FALSE)
    }
  }
  if (design=="2x2x4") {
    seqs  <- 2
    bkni  <- 1/4
    bk    <- 1
    dfe   <- parse(text="3*n-4", srcfile=NULL)
    if(robust) dfe <- parse(text="n-2", srcfile=NULL)
    dfRRe <- parse(text="n-2", srcfile=NULL)
    Emse  <- (s2wT + s2wR)/2
    cvec  <- c(1,1)
  }
  if (design=="2x2x3") { # TRT|RTR design or TRR|RTT
    seqs  <- 2
    bkni  <- 3/8
    bk    <- 1.5
    dfe   <- parse(text="2*n-3", srcfile=NULL)
    if(robust) dfe <- parse(text="n-2", srcfile=NULL)
    dfRRe <- parse(text="n/2-1", srcfile=NULL) # balanced design
    Emse  <- (s2wT + s2wR)/2
    cvec  <- c(2,1) # dummy
  }
  if (length(n)==1){
    # for unbalanced designs we divide the ns by ourself
    # to have smallest imbalance
    nv <- nvec(n=n, grps=seqs)
    if (nv[1]!=nv[length(nv)]){
      message("Unbalanced design. n(i)=", paste(nv, collapse="/"), " assumed.")
    }
    C2 <- sum(1/nv)*bkni
    n <- sum(nv)
  } else {
    # check length
    if (length(n)!=seqs) stop("n must be a vector of length=",seqs,"!")
    
    C2 <- sum(1/n)*bkni
    nv <- n
    n <- sum(n)
  }
  if (design=="2x2x3"){
    dfTT <- nv[1]-1
    dfRR <- nv[2]-1
    Emse <- (nv[1]*(2*s2wT+s2wR)/3+nv[2]*(s2wT+2*s2wR)/3)/n
  }

  # sd of the sample mean T-R (point estimator)
  sdm  <- sqrt(Emse*C2)
  mlog <- log(theta0)
  df   <- eval(dfe)
  if (design!="2x2x3"){
    dfRR <- eval(dfRRe)
    dfTT <- dfRR
  } 
  
  if(setseed) set.seed(123456)
  p <- .pwr.scABEL(mlog, sdm, C2, Emse, cvec, df, s2wR, dfRR, s2wT, dfTT, 
                   design, nsims, CVswitch, r_const, CVcap, 
                   ln_lBEL=log(theta1),ln_uBEL=log(theta2), alpha=alpha)
    
  as.numeric(p)
}

# --- working horse for power calculation for pure scABEL
.pwr.scABEL <- function(mlog, sdm, C2, Emse, cvec, df, 
                          s2wR, dfRR, s2wT, dfTT, design,
                          nsims, CVswitch=0.3, r_const=0.760, CVcap=0.5,
                          ln_lBEL=log(0.8), ln_uBEL=log(1.25), alpha=0.05)
{
  tcrit    <- qt(1-alpha,df)
  s2switch <- log(1.0 + CVswitch^2)
  s2cap    <- log(1.0 + CVcap^2)
  capABEL  <- sqrt(s2cap)*r_const
  
  denom    <- sum(cvec)
  
  # result vector
  counts        <- rep.int(0, times=4)
  names(counts) <- c("BE", "BEwl", "BEpe", "BEabe")
  # to avoid memory problems
  chunks <- 1
  nsi    <- nsims
  if (nsims>1E7) {
    chunks <- round(nsims/1E7,0)
    nsi    <- 1E7
  } 
  for (iter in 1:chunks){
    # simulate sample mean via its normal distribution
    means  <- rnorm(nsi, mean=mlog, sd=sdm)
    # simulate sample value s2wR via chi-square distri
    # shifting this after the sample mse sims will give
    # changes in the order of max. 4E-4 compared to V1.1-02!
    s2wRs  <- s2wR*rchisq(nsi, dfRR)/dfRR
    s2wTs <- s2wT*rchisq(nsi, dfTT)/(dfTT)
    if (design=="2x3x3"){
      mses  <- (cvec[1]*s2wTs + cvec[2]*s2wRs)/denom
    } 
    if (design=="2x2x4"){
      # simulate sample mse 
      mses  <- Emse*rchisq(nsi, df)/df
      #--- 'mean' of both attempts V1.1-02/V1.1-03 in case of 2x2x4
      mses  <- (mses + (cvec[1]*s2wTs + cvec[2]*s2wRs)/denom)/2
    }
    if (design=="2x2x3"){
      # simulate sample mse 
      mses  <- Emse*rchisq(nsi, df)/df
      # --- 'mean' of mses and mses calculated from components
      mses  <- 0.5*mses + 
               0.5*((dfTT+1)*(2*s2wTs + s2wRs)/3
                   +(dfRR+1)*(s2wTs + 2*s2wRs)/3)/(dfTT+dfRR+2)
      # from components
#      mses  <- ((dfTT+1)*(2*s2wTs + s2wRs)/3 +
#               (dfRR+1)*(s2wTs + 2*s2wRs)/3)/(dfTT+dfRR+2)
    }
    #--- EMA widened limits in log-domain
    uABEL   <- +sqrt(s2wRs)*r_const
    lABEL   <- -uABEL
    
    # no 'mixed' procedure, no cap
#     # use ABE limits if CV < CVswitch
#     lABEL[s2wRs<=s2switch] <- ln_lBEL
#     uABEL[s2wRs<=s2switch] <- ln_uBEL
#     # cap limits if CVcap not Inf
#     if (is.finite(CVcap)){
#       lABEL[s2wRs>s2cap] <- -capABEL
#       uABEL[s2wRs>s2cap] <- +capABEL
#     }
    #--- 90% CIs for T-R
    hw  <- tcrit*sqrt(C2*mses)
    lCL <- means - hw 
    uCL <- means + hw
    rm(hw)
    #--- 90% CI in 'widened' limits? 
    BE   <- (lABEL<=lCL) & (uCL<=uABEL)

    counts["BE"] <- counts["BE"]  + sum(BE)
  } # end over chunks
  
  # return the counts
  counts["BE"]/nsims
}