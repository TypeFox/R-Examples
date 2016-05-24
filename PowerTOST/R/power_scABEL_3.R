#---------------------------------------------------------------------------
# Simulate partial and full replicate design and scaled ABE power
# using 'widened' limits (EMA)
# 
# simulate via distributions of intra-subject contrasts (FDA) and afterwards
# estimation of mse 
# Author: dlabes
#---------------------------------------------------------------------------

# degrees of freedom for the TR/RR  analysis: 
# Using the intrasubject contrasts T-R and R-R and analyze them  
# by sequence groups the df's = n-seq (robust df's).
# 2x3x3  dfRR = n-3
# 2x2x4  dfRR = n-2
# 2x2x3  dfRR = n/2 - 2

power.scABEL3 <- function(alpha=0.05, theta1, theta2, theta0, CV, n,   
                          design=c("2x3x3", "2x2x4", "2x2x3"), 
                          regulator = c("EMA", "FDA"),
                          nsims=1E5, details=FALSE, setseed=TRUE)
{
  if (missing(CV)) stop("CV must be given!")
  if (missing(n))  stop("Number of subjects n must be given!")

  if (missing(theta0)) theta0 <- 0.95
  if (missing(theta1) & missing(theta2)) theta1 <- 0.8
  if (missing(theta2)) theta2 <- 1/theta1
  
  ptm <- proc.time()
  
  CVwT <- CV[1]
  if (length(CV)==2) CVwR <- CV[2] else CVwR <- CVwT
  s2wT <- log(1.0 + CVwT^2)
  s2wR <- log(1.0 + CVwR^2)

  regulator <- match.arg(regulator)
  CVswitch  <- 0.3
  r_const   <- 0.76 # or better log(theta2)/CV2se(0.3)? EMA & ANVISA
  if (regulator=="FDA") r_const <- log(1.25)/0.25 # or better log(theta2)/0.25?
  CVcap <- 0.5
  if (regulator=="FDA") CVcap=Inf
  
  design <- match.arg(design)
  if (design=="2x3x3") {
    seqs <- 3
    bkni <- 1/6
    bk   <- 1.5
    dfe   <- parse(text="n-3", srcfile=NULL)
    dfCIe <- parse(text="2*n-3", srcfile=NULL)
    dfRRe <- parse(text="n-3", srcfile=NULL)
    # according to McNally et al.
    # verified via simulations:
    Es2I  <- s2wT + s2wR/2
    Emse  <- (s2wT + 2.0*s2wR)/3
    # warning in case of CVwR != CVwT
    if (abs(s2wT-s2wR)>1e-4){
      warning(paste("Heteroscedasticity in the 2x3x3 design may led", 
                    "to poor accuracy of power!"), call.=FALSE)
    }
  }
  if (design=="2x2x4") {
    seqs  <- 2
    bkni  <- 1/4
    bk    <- 1
    dfe   <- parse(text="n-2", srcfile=NULL)
    dfCIe <- parse(text="3*n-4", srcfile=NULL)
    dfRRe <- parse(text="n-2", srcfile=NULL)
    # sd^2 of the differences T-R from their components
    Es2I  <- (s2wT + s2wR)/2 
    Emse  <- (s2wT + s2wR)/2 
  }
  if (design=="2x2x3") {
    seqs  <- 2
    bkni  <- 3/8
    bk    <- 1.5
    dfe   <- parse(text="n-2", srcfile=NULL)
    dfCIe <- parse(text="2*n-3", srcfile=NULL)
    # next was pre-V1.2-08
#     dfRRe <- parse(text="n/2-2", srcfile=NULL) # for balanced designs
#     dfTTe <- parse(text="n/2-2", srcfile=NULL) # for balanced designs
    # correct should be (only 1 sequence for each, f.i. RR from RTR):
    dfRRe <- parse(text="n/2-1", srcfile=NULL) # for balanced designs
    dfTTe <- parse(text="n/2-1", srcfile=NULL) # for balanced designs
    # sd^2 of the differences T-R from their components
    Es2I  <- 1.5*(s2wT + s2wR)/2               # for balanced design 
    Emse  <- (s2wT + s2wR)/2
  }
  
  if (length(n)==1){
    # then we assume n=ntotal
    # for unbalanced designs we divide the ns by ourself
    # to have only small imbalance
    nv <- nvec(n=n, grps=seqs)
    if (nv[1]!=nv[length(nv)]) {
      message("Unbalanced design. n(i)=", paste(nv, collapse="/"), " assumed.")
    } 
  } else {
    # then we assume n = vector of n's in sequences
    # check length
    if (length(n)!=seqs) stop("n must be a vector of length=",seqs,"!")
    nv <- n
  }
  C2 <- sum(1/nv)*bkni
  n  <- sum(nv)
  
  df   <- eval(dfe)
  dfCI <- eval(dfCIe)

  if (design=="2x2x3"){
    dfTT <- nv[1]-1
    dfRR <- nv[2]-1
    # where does this next came from? McNally?
    Es2I <- (dfRR*(s2wT + s2wR/2)+dfTT*(s2wT/2 + s2wR))/(dfRR+dfTT)
    Emse <- (nv[1]*(2*s2wT+s2wR)/3+nv[2]*(s2wT+2*s2wR)/3)/n
    # here we leave the algorithm as was
    Es2I <- Emse
    # warning in case of unbalanced design and heteroscdasticity
    # if (abs(s2wT - s2wR)>1e-5 & abs(dfRR-dfTT)>2){
    #   warning(paste("Heteroscedasticity in unbalanced 2x2x3 design may led",
    #           "to poor accuracy of power!"), call.=FALSE)
    # }
  } else {
    dfRR <- eval(dfRRe)
    dfTT <- dfRR
  }
  #cat("dfRR=", dfRR," dfTT=",dfTT," E(s2I)=", Emse, "\n")
  # sd of the mean T-R (point estimator)
  sdm  <- sqrt(Emse*C2)
  mlog <- log(theta0)
  
  if(setseed) set.seed(123456)
  p <- .pwr.scABEL.ISC2(mlog, sdm, C2, Es2I, df, dfCI, s2wR, dfRR, s2wT, dfTT,
                        design=design, nsims=nsims, 
                        CVswitch=CVswitch, r_const=r_const, CVcap=CVcap, 
                        ln_lBEL=log(theta1),ln_uBEL=log(theta2), alpha=alpha)
    
  if (details) {
    ptm <- summary(proc.time()-ptm)
    message(nsims,"sims. Time elapsed (sec): ", 
            formatC(ptm["elapsed"], digits=2), "\n")
    #print(ptm)
    # return also the components
    names(p) <- c("p(BE)", "p(BE-wABEL)", "p(BE-pe)", "p(BE-ABE)")
    p
  } else {
    # return only the 'power'
    as.numeric(p["BE"])
  }
}

# working horse 
.pwr.scABEL.ISC2 <- function(mlog, sdm, C2, Es2I, df, dfCI, s2wR, dfRR, s2wT, dfTT,
                             design, nsims, CVswitch=0.3, r_const=0.892574, CVcap=0.5, 
                             ln_lBEL=log(0.8), ln_uBEL=log(1.25), alpha=0.05)
{
  #tval     <- qt(1-alpha,df)
  tval     <- qt(1-alpha,dfCI)    # usual df's for ANOVA for CI
  s2switch <- log(CVswitch^2+1)
  s2cap    <- log(CVcap^2+1)
  
  counts <- rep.int(0, times=4)
  names(counts) <- c("BE", "BEwl", "BEpe", "BEabe")
  # to avoid memory problems for high number of sims
  chunks <- 1
  nsi    <- nsims
  if (nsims>1E7) {
    chunks <- round(nsims/1E7,0)
    nsi    <- 1E7
  } 
  for (iter in 1:chunks) {
    # same order of sims to avoid samll differences compared to previous version
    # simulate sample mean via its normal distribution
    means  <- rnorm(nsi, mean=mlog, sd=sdm)
    # simulate sample value s2wRs via chi-square distri
    s2wRs  <- s2wR*rchisq(nsi, dfRR)/dfRR
    # simulate sample value s2wTs via chi-square distri
    s2wTs  <- s2wT*rchisq(nsi, dfTT)/dfTT
    
    if (design=="2x3x3") {
      # simulate sample s2Is via chi-square distri
      s2Is   <- Es2I*rchisq(nsi, df)/df
      #we have s2wT + s2wR/2 and need (s2wT+2*s2wR)/3
      mses <- s2Is/3 + 0.5*s2wRs
    } 
    if (design=="2x2x4") {
      # simulate sample s2Is via chi-square distri
      s2Is   <- Es2I*rchisq(nsi, df)/df
      # we make a second estimate of mse and take the mean
      # thus we may achieve the ANOVA df's, at least approximately 
      mses <- (s2Is + 0.5*(s2wTs + s2wRs))/2
    }  
    if (design=="2x2x3") {
      # at moment can't figure out how to do it here thus leave the algo as was
      mses  <- Es2I*rchisq(nsi, dfCI)/dfCI
      # dfTT+1 = n1, dfRR = n2
      mses <- 0.5*mses +
              0.5*((dfTT+1)*(2*s2wTs + s2wRs)/3
                  +(dfRR+1)*(s2wTs + 2*s2wRs)/3)/(dfTT+dfRR+2)
    }
    SEs <- sqrt(mses*C2)
    # conventional (1-2*alpha) CI's for T-R
    hw  <- tval*SEs
    lCL <- means - hw 
    uCL <- means + hw
    # conventional ABE
    BEABE <- (lCL>=ln_lBEL) & (uCL<=ln_uBEL)
    
    #--- widened limits in log-domain
    uABEL <- +sqrt(s2wRs)*r_const
    # cap on 'widened' limits
    uABEL[s2wRs>s2cap] <- sqrt(s2cap)*r_const
    # BE using widened acceptance limits
    BE <- (lCL>=-uABEL) & (uCL<=uABEL)
    
    # pe constraint
    BEpe  <- ( means>=ln_lBEL & means<=ln_uBEL )
    
    # save memory
    rm(SEs, hw, uABEL)
    
    # if CV < CV switch use ABE, else scABEL
    BE    <- ifelse(s2wRs>s2switch, BE, BEABE)

    counts["BEabe"] <- counts["BEabe"] + sum(BEABE)
    counts["BEpe"]  <- counts["BEpe"]  + sum(BEpe)
    counts["BEwl"]  <- counts["BEwl"]  + sum(BE)
    counts["BE"]    <- counts["BE"]    + sum(BE & BEpe) # with pe constraint
    
  } # end over chunks
  
  # return the pBEs
  return(counts/nsims)
}
