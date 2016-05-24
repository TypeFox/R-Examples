# --------------------------------------------------------------------------
# power (or alpha) of 2-stage studies according to Potvin et. al. 
# methods "B" and "C", modified to include a futility criterion Nmax,
# modified to use PE and mse of stage 1 in power calculation steps as well 
# as in sample size estimation, Karalis/Macheras modifications
#
# author D.L.
# --------------------------------------------------------------------------
# require(PowerTOST)

power.2stage.KM <- function(method=c("C","B"), alpha0=0.05, alpha=c(0.0294,0.0294),
                            n1, CV, targetpower=0.8, pmethod=c("nct","exact"),
                            Nmax=150, theta0, theta1, theta2, 
                            npct=c(0.05, 0.5, 0.95), nsims, setseed=TRUE, 
                            details=FALSE)
{
  if (missing(CV)) stop("CV must be given!")
  if (CV<=0)       stop("CV must be >0!")
  
  if (missing(n1)) stop("Number of subjects in stage 1 must be given!")
  if (n1<=0)       stop("Number of subjects in stage 1 must be >0!")
  
  if (missing(theta1) & missing(theta2))  theta1 <- 0.8
  if (!missing(theta1) & missing(theta2)) theta2 <- 1/theta1
  if (missing(theta1) & !missing(theta2)) theta1 <- 1/theta2
  
  if (missing(theta0)) theta0 <- 0.95
  
  if (n1>Nmax) stop("n1>Nmax doestn't make sense!")
  
  if(missing(nsims)){
    nsims <- 1E5
    if(theta0<=theta1 | theta0>=theta2) nsims <- 1E6
  }
  
  # check if Potvin B or C
  method  <- match.arg(method)
  # check if power calculation method is nct or exact
  pmethod <- match.arg(pmethod)
  
  if(details){
    cat(nsims,"sims. Stage 1")
  }
  # start timer
  ptm  <- proc.time()
  
  if (setseed) set.seed(1234567)

  ltheta1 <- log(theta1)
  ltheta2 <- log(theta2)
  mlog    <- log(theta0)
  mse     <- CV2mse(CV)
  bk      <- 2   # 2x2x2 crossover design const
  
  # reserve memory
  BE    <- rep.int(NA, times=nsims)
# ----- stage 1 ----------------------------------------------------------
  Cfact <- bk/n1
  df    <- n1-2
  tval  <- qt(1-alpha[1], df)
  sdm   <- sqrt(mse*Cfact)
  # simulate point est. via normal distribution
  pes   <- rnorm(n=nsims, mean=mlog, sd=sdm)
  # simulate mse via chi-squared distribution
  mses  <- mse*rchisq(n=nsims, df=df)/df
  
  # K&M have the test pe in BE range first
  # may speed up somewhat
  # next construction defines some zone at the BE acceptance limits
  # where sample size is practical Inf, these are counted as outside
  outside <- ((pes-ltheta1)<1.25e-5 | (ltheta2-pes)<1.25e-5)
  BE      <- !outside  # =FALSE for outside -> FAIL
  BE[BE==TRUE] <- NA   # not outside, not yet decided
  
  if(method=="C"){
    mses_tmp <- mses[is.na(BE)]
    pes_tmp  <- pes[is.na(BE)]
    # if method=C then calculate power for alpha0=0.05, mse and pe from stage 1
    pwr <- mapply(.calc.power, diffm=pes_tmp, sem=sqrt(bk*mses_tmp/n1),
                  MoreArgs=list(alpha=alpha0, ltheta1=ltheta1, ltheta2=ltheta2,
                                df=df, method=pmethod))
    tval0 <- qt(1-alpha0, df)
    hw    <- tval0*sqrt(Cfact*mses_tmp)
    lower <- pes_tmp - hw
    upper <- pes_tmp + hw
    # fail or pass
    BE0    <- lower>=ltheta1 & upper<=ltheta2
    # if power>0.8 then calculate CI for alpha=0.05
    # i.e. if power<0.8 then 
    BE0[pwr<targetpower] <- NA # not yet decided
    # combine these with the previous
    BE[is.na(BE)] <- BE0
    rm(BE0)
  }
  # method "B" or power<=0.8 in method "C":
  # evaluate BE (CI) with alpha=alpha1
  mses_tmp <- mses[is.na(BE)]
  pes_tmp  <- pes[is.na(BE)]
  BE1 <- rep.int(NA, times=length(mses_tmp))
  hw    <- tval*sqrt(Cfact*mses_tmp)
  lower <- pes_tmp - hw
  upper <- pes_tmp + hw
  rm(hw)
  BE1   <- lower>=ltheta1 & upper<=ltheta2
  # take care of memory
  rm(lower, upper)
  if (method=="C"){
    #if BE met -> PASS stop
    #if not BE -> goto sample size estimation i.e flag BE1 as NA
    BE1[!BE1] <- NA
  } else { 
    # method B
    # evaluate power at alpha[1] using PE and mse of stage 1
    pwr <- mapply(.calc.power, diffm = pes_tmp, sem = sqrt(bk*mses_tmp/n1), 
                  MoreArgs = list(alpha = alpha[1], ltheta1 = ltheta1, 
                                  ltheta2 = ltheta2, df = df, method = pmethod))
    # if BE met then decide BE regardless of power
    # if not BE and power<0.8 then goto stage 2
    BE1[ !BE1 & pwr<targetpower ] <- NA 
    # take care of memory
    rm(pwr)
  }
  # combine 'stage 0' from method C and stage 1
  BE[is.na(BE)] <- BE1
  # take care of memory, done with stage 1
  rm(BE1)
  
  if(details){
    cat(" - Time consumed (secs):\n")
    print(round((proc.time()-ptm),1))
  }

  # ------sample size for stage 2 -----------------------------------------
  ntot     <- rep(n1, times=nsims)
  stage    <- rep(1, times=nsims)
  # filter out those were stage 2 is necessary
  pes_tmp  <- pes[is.na(BE)]
  
  # Maybe we are already done with stage 1
  if(length(pes_tmp)>0){
    if(details){
      cat("Keep calm. Sample sizes (", length(pes_tmp),
          " studies) for stage 2\n", sep="")
      cat("will be estimated. May need some time.\n")
    }
    # preliminary setting stage=2 for those not yet decided BE
    # may be altered for those with nt>Nmax or nt=Inf 
    # from sample size est. if pe outside acceptance range
    # see below
    stage[is.na(BE)] <- 2
    mses_tmp <- mses[is.na(BE)]
    BE2      <- rep.int(NA, times=length(mses_tmp))
    s2       <- rep.int(2, times=length(mses_tmp))
    #------ sample size for stage 2 ---------------------------------------
    ptms <- proc.time()
    # use mse1 & pe1 as described in Karalis/Macheras
    # sample size function returns Inf if pe1 is outside acceptance range
#     nt <- mapply(FUN=.sampleN, mse=mses_tmp, ltheta0=pes_tmp, 
#                  MoreArgs=list(alpha=alpha[2], targetpower=targetpower, 
#                                ltheta1=ltheta1, ltheta2=ltheta2, 
#                                method=pmethod, bk=bk))
    nt <- .sampleN2(alpha=alpha[2], targetpower=targetpower, ltheta0=pes_tmp,
                    mse=mses_tmp, ltheta1=ltheta1, ltheta2=ltheta2, 
                    method=pmethod, bk=bk)
    n2  <- ifelse(nt>n1, nt - n1, 0)
    
    if(details){
      if(nsims<=1E5 & pmethod!="exact"){
        cat("Time consumed (secs):\n")
        print(round((proc.time()-ptms),1))
      } else {
        cat("Time consumed (min):\n")
        print(round((proc.time()-ptms)/60,2))
      }
    }
    # futility rule: if nt > Nmax -> stay with stage 1 result not BE
    # ntotal = n1 reasonable?
    if (is.finite(Nmax) | any(!is.finite(nt))){
      # sample size may return Inf if PE is used in ss estimation
      # in that case we stay with stage 1 result
      BE2[!is.finite(n2) | (n1+n2)>Nmax] <- FALSE
      # and we are counting these for stage 1
      s2[BE2==FALSE]  <- 1
      # debug print
      # cat(sum(!BE2, na.rm=T)," cases with nt>Nmax or nt=Inf\n")
      # save 
      stage[is.na(BE)] <- s2
      # save the FALSE and NA in BE
      BE[is.na(BE)]    <- BE2
      # filter out those were BE was yet not decided
      pes_tmp  <- pes_tmp[is.na(BE2)]
      mses_tmp <- mses_tmp[is.na(BE2)]
      n2       <- n2[is.na(BE2)]
    }
    # ---------- stage 2 evaluation --------------------------------------
    m1    <- pes_tmp
    SS1   <- (n1-2)*mses_tmp
    nsim2 <- length(pes_tmp)
    # to avoid warnings for n2=0 in rnorm() and rchisq()
    ow    <- options("warn")
    options(warn=-1)
    m2    <- ifelse(n2>0, rnorm(n=nsim2, mean=mlog, sd=sqrt(mse*bk/n2)), 0)
    SS2   <- ifelse(n2>2, (n2-2)*mse*rchisq(n=nsim2, df=n2-2)/(n2-2), 0)
    # reset options
    options(ow) 
    SSmean <- ifelse(n2>0, (m1-m2)^2/(2/n1+2/n2), 0)
    nt     <- n1+n2
    df2    <- ifelse(n2>0, nt-3, n1-2)
    pe2    <- ifelse(n2>0, (n1*m1+n2*m2)/nt, pes_tmp)
    mse2   <- ifelse(n2>0, (SS1+SSmean+SS2)/df2, mses_tmp)
    # take care of memory
    rm(m1, m2, SS1, SS2, SSmean)
    # calculate CI for stage 2 with alpha[2]
    tval2 <- qt(1-alpha[2], df2)
    hw    <- tval2*sqrt(mse2*bk/nt)
    lower <- pe2 - hw
    upper <- pe2 + hw
    BE2   <- lower>=ltheta1 & upper<=ltheta2
    # combine stage 1 & stage 2 BE statements
    ntot[is.na(BE)]  <- nt
    BE[is.na(BE)]    <- BE2
    # done with them
    rm(BE2, nt, lower, upper, hw)
  } # end stage 2 calculations
  # take care of memory - may be superflous here since function ends
  rm(pes_tmp, mses_tmp)
  # the return list
  res <- list(design="2x2 crossover", method=method, modified="KM",
              alpha0=ifelse(method=="C",alpha0,NA), alpha=alpha, 
              CV=CV, n1=n1, targetpower=targetpower, pmethod=pmethod, 
              theta0=exp(mlog), theta1=theta1, theta2=theta2, Nmax=Nmax, 
              nsims=nsims,
              # results 
              pBE=sum(BE)/nsims, pBE_s1=sum(BE[stage==1])/nsims, 
              pct_s2=100*length(BE[stage==2])/nsims, 
              nmean=mean(ntot), nrange=range(ntot), nperc=quantile(ntot, p=npct)
              )

  # table object summarizing the discrete distri of ntot
  # only given back if usePE=FALSE or if usePE=TRUE then Nmax must be finite
  # since here usePE not available but is used as TRUE this reduces to is.finite(Nmax)
  if (is.finite(Nmax)){
    res$ntable <- table(ntot)
  }
  
  # output is now done via S3 print method

  class(res) <- c("pwrtsd", "list")
  return(res)
  
} #end function
