# --------------------------------------------------------------------------
# power (or alpha) of 2-stage group sequential 2x2 crossover studies 
# with no sample size adaption (n1, n2 predefined)
#
# author D.Labes
# --------------------------------------------------------------------------
# require(PowerTOST)

power.2stage.GS <- function(alpha=c(0.0294,0.0294), n, CV, theta0, theta1, 
                            theta2,  fCrit=c("PE","CI"), fClower, fCupper, 
                            nsims, setseed=TRUE, details=FALSE)
{
  if (missing(CV)) stop("CV must be given!")
  if (CV<=0)       stop("CV must be >0!")
  
  if (missing(n))   stop("Number of subjects in the two stages must be given.")
  if (length(n)!=2) stop("n must be a vector of length 2.")
  if (any(n<=2))    stop("Number of subjects in stages must be >2.")
  
  if (missing(theta1) & missing(theta2))  theta1 <- 0.8
  if (!missing(theta1) & missing(theta2)) theta2 <- 1/theta1
  if (missing(theta1) & !missing(theta2)) theta1 <- 1/theta2
  
  if (missing(theta0)) theta0 <- 0.95

  if(missing(nsims)){
    if(theta0<=theta1 | theta0>=theta2) nsims <- 1E6 else  nsims <- 1E5
  }
  
  fCrit <- match.arg(fCrit)
  if (missing(fClower) & missing(fCupper))  fClower <- 0 # no futility crit.
  if (fClower<0) fClower <- 0 # silently correct it
  if (missing(fClower) & !missing(fCupper)) fClower <- 1/fCupper
  if (!missing(fClower) & missing(fCupper)) fCupper <- 1/fClower

  
  if(details){
    cat(nsims, "sims. Stage 1")
  }
  # start timer
  ptm  <- proc.time()
  
  if (setseed) set.seed(1234567)

  ltheta1 <- log(theta1)
  ltheta2 <- log(theta2)
  mlog    <- log(theta0)
  mse     <- CV2mse(CV)
  bk      <- 2   # 2x2x2 crossover design const
  lfC1   <- log(fClower)
  lfC2   <- log(fCupper)
  # reserve memory
  BE     <- rep.int(TRUE, times=nsims)
  
# ----- stage 1 ----------------------------------------------------------
  Cfact <- bk/n[1]
  df    <- n[1]-2
  tval  <- qt(1-alpha[1], df)
  sdm   <- sqrt(mse*Cfact)
  # simulate point est. via normal distribution
  pes   <- rnorm(n=nsims, mean=mlog, sd=sdm)
  # simulate mse via chi-squared distribution
  mses  <- mse*rchisq(n=nsims, df=df)/df
  
  # calculate 1-2*alpha CIs
  hw <- tval*sqrt(Cfact*mses)
  lower <- pes - hw
  upper <- pes + hw
  # fail or pass
  BE   <- lower>=ltheta1 & upper<=ltheta2
  rm(hw)

  # ----- stage 2 ---------------------------------------------------------
  # only those not BE in stage 1
  stage <- rep.int(1, times=length(BE))
  stage[BE==FALSE] <- 2
  # any power criterion? aka futility criterion?
  if (fCrit=="PE"){
    # rule out those where point est. is not in acceptance range
    # or not in a prespecified range, f.i. 0.85 ... 1/0.85=1.176471
    s2 <- (pes>=lfC1 & pes<=lfC2) & stage==2
  } else {
    # another possibility is Gould (phi=1):
    # rule out those with CI totally outside acceptance range
    s2 <- !(lower>lfC2 | upper<lfC1) & stage==2
    # or not in a prespecified range, f.i. 0.85 ... 1/0.85=1.176471
    # Gould phi=0 reads: CI doesnt contain zero (1 in orginal domain)
  }
  stage[s2==TRUE]  <- 2
  stage[s2==FALSE] <- 1
  
  ns <- rep.int(1, times=length(BE))
  ns[s2==FALSE] <- n[1]
  ns[s2==TRUE]  <- n[1]+n[2]
  
  nsims2  <- sum(s2)

  # time for stage 1
  if(details){
    cat(" - Time consumed (sec):\n")
    print(proc.time()-ptm)
  }

  if (nsims2>0){
    m1    <- pes[s2==TRUE]
    SS1   <- (n[1]-2)*mses[s2==TRUE]
    # --- simulate stage 2 data
    Cfact <- bk/n[2]
    sdm   <- sqrt(mse*Cfact)
    # simulate point est. via normal distribution
    m2    <- rnorm(n=nsims2, mean=mlog, sd=sdm)
    # simulate mse/SS2 via chi-squared distribution
    SS2   <- 0
    if(n[2]>2) SS2 <- mse*rchisq(n=nsims2, df=n[2]-2)
    # --- poole over stage 1 & stage 2 data
    # pool according to Potvin et al.
    ntot <- sum(n)
    df2  <- ntot-3 # including stage
    tval2 <- qt(1-alpha[2], df2)
    # The Canada guidance talks: "This method precludes the need for 
    # a stage effect in the model" -> setting SSmean to zero?
    SSmean <- (m1-m2)^2/(2/n[1]+2/n[2])
    #SSmean <- 0
    # mean stage 1 + stage 2
    pe2 <- (n[1]*m1+n[2]*m2)/ntot
    # mse stage 1 + stage 2 with stage as covariate
    mse2   <- (SS1+SSmean+SS2)/df2
    # take care of memory
    rm(m1, m2, SS1, SS2, SSmean)
    # calculate 1-2*alpha CI for stage 2 with alpha[2]
    hw    <- tval2*sqrt(mse2*bk/ntot)
    lower <- pe2 - hw
    upper <- pe2 + hw
    BE2   <- lower>=ltheta1 & upper<=ltheta2
    # take care of memory
    rm(hw, lower, upper)
    # --- combine BE and BE2
    BE[s2==TRUE] <- BE2
  }

  # the return list
  res <- list(alpha=alpha, CV=CV, n=n, theta0=exp(mlog), theta1=theta1, 
              theta2=theta2, fCrit=fCrit, fCrange=c(fClower, fCupper),
              nsims=nsims,
              # results 
              pBE=sum(BE)/nsims, pBE_s1=sum(BE[stage==1])/nsims,
              pct_s2=100*length(BE[stage==2])/nsims
              # mean N?, median?
              )
  
  # histogram data or table data are not available here
  if (details){
    cat("Total time consumed (sec):\n")
    print(proc.time()-ptm)
    cat("\n")
  }
  
  # output is now via S3 print method

  class(res) <- c("pwrtsd", "list")
  return(res)
  
} #end function
