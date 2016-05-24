# --------------------------------------------------------------------------
# power (or alpha) of 2-stage studies according to Potvin et. al. 
# method "B" or "C" modified to include a futility criterion for PE or CI and
# modified to use PE of stage 1 in sample size estimation
#
# author D.L.
# --------------------------------------------------------------------------
# require(PowerTOST)
# source("./R/sampsiz2.R")
# source("./R/sampsiz_n0.R")
# source("./R/power.R")

power.2stage.fC <- function(method=c("B", "C"), alpha0=0.05, alpha=c(0.0294,0.0294), 
                            n1, CV, GMR, targetpower=0.8, 
                            pmethod=c("nct","exact", "shifted"), 
                            usePE=FALSE, powerstep=TRUE, min.n2=0, max.n=Inf,
                            fCrit=c("CI", "PE"), fClower, fCupper, theta0, theta1, 
                            theta2, npct=c(0.05, 0.5, 0.95), nsims, setseed=TRUE, 
                            details=FALSE)
{ 
  # chek method, check futility criterion
  method <- match.arg(method)
  fCrit  <- match.arg(fCrit)
  
  if (missing(CV)) stop("CV must be given!")
  if (CV<=0)       stop("CV must be >0!")
  
  if (missing(n1)) stop("Number of subjects in stage 1 (n1) must be given!")
  if (n1<=0)       stop("Number of subjects in stage 1 (n1) must be >0!")
  
  if (missing(GMR)) GMR <- 0.95
  
  if (missing(theta1)  & missing(theta2))  theta1 <- 0.8
  if (!missing(theta1) & missing(theta2))  theta2 <- 1/theta1
  if (missing(theta1)  & !missing(theta2)) theta1 <- 1/theta2
  
  if (GMR<=theta1 | GMR>=theta2) stop("GMR must be within acceptance range!")
  
  if (missing(theta0)) theta0 <- GMR
  
  if(missing(nsims)){
    if(theta0<=theta1 | theta0>=theta2) nsims <- 1E6 else  nsims <- 1E5
  }
  
  if(fCrit=="PE"){
    if (missing(fClower) & missing(fCupper))  fClower <- 0.8
    if (missing(fClower) & !missing(fCupper)) fClower <- 1/fCupper
    if (!missing(fClower) & missing(fCupper)) fCupper <- 1/fClower
  }
  if(fCrit=="CI"){
    if (missing(fClower) & missing(fCupper))  fClower <- 0.925
    if (missing(fClower) & !missing(fCupper)) fClower <- 1/fCupper
    if (!missing(fClower) & missing(fCupper)) fCupper <- 1/fClower
  }

  if(min.n2!=0 & min.n2<2) stop("min.n2 has to be at least +2 if >0.")
  # make even (round up)
  if( min.n2%%2 != 0) {
    min.n2 <- min.n2 + min.n2%%2
    message("min.n2 rounded to even", min.n2)
  }
  
  # check max.n
  if (n1>=max.n) stop("max.n <= n1 doestn't make sense!")
  
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
  ln_fClower <- log(fClower)
  ln_fCupper <- log(fCupper)
  lGMR    <- log(GMR)
  mlog    <- log(theta0)
  mse     <- CV2mse(CV)
  bk      <- 2   # 2x2x2 crossover design const
  # reserve memory
  BE      <- rep.int(NA, times=nsims)
  
# ----- stage 1 ----------------------------------------------------------
  Cfact <- bk/n1
  df    <- n1-2
  tval  <- qt(1-alpha[1], df)
  tval0 <- qt(1-alpha0, df)

  sdm   <- sqrt(mse*Cfact)
  # simulate point est. via normal distribution
  pes   <- rnorm(n=nsims, mean=mlog, sd=sdm)
  # simulate mse via chi-squared distribution
  mses  <- mse*rchisq(n=nsims, df=df)/df
  if(method=="C"){
    # if method=C then calculate power for alpha0=0.05 and plan GMR
    pwr <- .calc.power(alpha=alpha0, ltheta1=ltheta1, ltheta2=ltheta2, 
                       diffm=lGMR, sem=sqrt(bk*mses/n1), df=df,  
                       method=pmethod)
    
    tval0 <- qt(1-alpha0, df)
    hw    <- tval0*sqrt(Cfact*mses)
    lower <- pes - hw
    upper <- pes + hw
    # fail or pass
    BE    <- lower>=ltheta1 & upper<=ltheta2
    # if power>0.8 then calculate CI for alpha=0.05
    # i.e. if power<0.8 then 
    BE[pwr<targetpower] <- NA # not yet decided
    powerstep <- TRUE
    # take care of memory
    rm(hw, lower, upper, pwr)
  }
  mses_tmp <- mses[is.na(BE)]
  pes_tmp  <- pes[is.na(BE)]
  BE1 <- rep.int(NA, times=length(mses_tmp))
  # calculate 1-2*alpha CI for alpha=alpha1
  hw    <- tval*sqrt(Cfact*mses_tmp)
  lower <- pes_tmp - hw
  upper <- pes_tmp + hw
  BE1   <- lower>=ltheta1 & upper<=ltheta2
  # browser()
  if (method=="C"){
    #if BE met -> PASS stop
    #if not BE -> goto sample size estimation i.e flag BE1 as NA
    BE1[!BE1] <- NA
  } else { 
    # method B / E
    if(powerstep){
      # evaluate power at alpha[2] according to Xu et al. Method E, F
      pwr <- .calc.power(alpha=alpha[2], ltheta1=ltheta1, ltheta2=ltheta2, 
                         diffm=lGMR, sem=sqrt(bk*mses_tmp/n1), df=df, 
                         method=pmethod)
      # Potvin method E:
      # if not BE and if power >= 0.8 (targetpower) make a second BE evaluation 
      # with alpha[2]
      # but only if alpha[1] != alpha[2] ?
      BE12  <- BE1 # reserve memory
      BE11  <- BE1
      # BE decision at alpha[2]
      tval  <- qt(1-alpha[2], df)
      hw    <- tval*sqrt(Cfact*mses_tmp)
      lower <- pes_tmp - hw
      upper <- pes_tmp + hw
      BE12  <- lower>=ltheta1 & upper<=ltheta2
      # browser()
      # if BE(a1) then BE1=TRUE, regardless of power
      BE1[BE11==TRUE] <- TRUE 
      # if not BE(a1) but power >= 0.8 then make BE decision at alpha2
      BE1[BE11==FALSE & pwr>=targetpower] <- BE12[BE11==FALSE & pwr>=targetpower]
      # if not BE(a1) and power<0.8 then not decided (marker NA)
      # will be further decided by futility criterion
      BE1[BE11==FALSE & pwr<targetpower]  <- NA 
      # keep care of memory
      rm(BE11, BE12)
    } else {
      # we do not calculate power
      BE1[ !BE1 ] <- NA  # not decided yet -> stage 2
    }
  }
  # ----- check pe /ci outside futility range
  if(fCrit=="PE"){
    outside <- ((pes_tmp-ln_fClower)<1.25e-5 | (ln_fCupper-pes_tmp)<1.25e-5)
  } else {
    # 90% (!) CI outside
    tval0   <- qt(1-alpha0, df)
    hw      <- tval0*sqrt(Cfact*mses_tmp)
    lower   <- pes_tmp - hw
    upper   <- pes_tmp + hw
    outside <- (lower > ln_fCupper) | (upper<ln_fClower)
  }
  BE1[is.na(BE1) & outside] <- FALSE
  # combine BE and BE1
  BE[is.na(BE)] <- BE1
  # take care of memory, done with them
  rm(hw, lower, upper, BE1, outside)
  
  # time for stage 1
  if(details){
    cat(" - Time consumed (secs):\n")
    print(round((proc.time()-ptm),1))
  }

  # ------ sample size for stage 2 -----------------------------------------
  ntot     <- rep(n1, times=nsims)
  stage    <- rep(1, times=nsims)
  # filter out those were stage 2 is necessary
  pes_tmp  <- pes[is.na(BE)]
  mses_tmp <- mses[is.na(BE)]
  
  # Maybe we are already done with stage 1
  if (length(pes_tmp)>0) {
    if (details) {
      cat("Keep calm. Sample sizes for stage 2 (", length(pes_tmp),
          " studies)\n", sep="")
      cat("will be estimated. May need some time.\n")
    }
    # preliminary setting stage=2 for those not yet decided BE
    # may be altered for those with nt>Nmax or nt=Inf 
    # from sample size est. if pe outside acceptance range
    # see below
    stage[is.na(BE)] <- 2
    BE2      <- rep.int(NA, times=length(mses_tmp))
    s2       <- rep.int(2, times=length(mses_tmp))
    #------ sample size for stage 2 ---------------------------------------
    ptms <- proc.time()
    if (usePE){
      # use mse1 & pe1 like in the paper of Karalis/Macheras
      # sample size function returns Inf if pe1 is outside acceptance range
#       nt <- mapply(FUN=.sampleN, mse=mses_tmp, ltheta0=pes_tmp, 
#                    MoreArgs=list(alpha=alpha[2], targetpower=targetpower, 
#                                  ltheta1=ltheta1, ltheta2=ltheta2,
#                                  method=pmethod, bk=2))
      nt <- .sampleN2(alpha=alpha[2], targetpower=targetpower, ltheta0=pes_tmp,
                      mse=mses_tmp, ltheta1=ltheta1, ltheta2=ltheta2, 
                      method=pmethod)
    } else {
      # use mse1 & GMR to calculate sample size (original Potvin)
#       nt <- mapply(FUN=.sampleN, mse=mses_tmp, 
#                    MoreArgs=list(alpha=alpha[2], targetpower=targetpower, 
#                                  ltheta0=lGMR, ltheta1=ltheta1, ltheta2=ltheta2,
#                                  method=pmethod, bk=2))
      nt <- .sampleN2(alpha=alpha[2], targetpower=targetpower, ltheta0=lGMR,
                      mse=mses_tmp, ltheta1=ltheta1, ltheta2=ltheta2, 
                      method=pmethod)
    }
    # cap to max.n
    nt[nt>max.n] <- max.n
    # n2
    n2 <- ifelse(nt>n1, nt - n1, 0)
    # assure a min.n2
    n2[n2<min.n2] <- min.n2 
    
    
    if(details){
      cat("Time consumed (secs):\n")
      print(round((proc.time()-ptms),1))
    }
    
    if (any(!is.finite(nt))){
      # sample size may return Inf if PE is used in ss estimation
      # for cases where PE is outside theta1 ... theta2
      # in that case we stay with stage 1
      BE2[!is.finite(n2)] <- FALSE
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
    # now simulate PE2 and SS2
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
    # combine stage 1 & stage 2
    ntot[is.na(BE)] <- nt
    BE[is.na(BE)]   <- BE2
    # done with them
    rm(BE2, nt, lower, upper, hw)
  } # ----- end of stage 2 calculations ----------------------------------
  # take care of memory
  rm(pes_tmp, mses_tmp)
  # the return list
  res <- list(design="2x2 crossover",
              method=method, alpha0=alpha0, # alpha0 is in method B also used for CI futility
              alpha=alpha, CV=CV, n1=n1, GMR=GMR, 
              targetpower=targetpower, pmethod=pmethod, 
              theta0=exp(mlog), theta1=theta1, theta2=theta2, usePE=usePE, 
              powerstep=powerstep, min.n2=min.n2, max.n=max.n,
              fCrit=fCrit, fCrange=c(fClower, fCupper), 
              nsims=nsims,
              # results 
              pBE=sum(BE)/nsims, pBE_s1=sum(BE[ntot==n1])/nsims,
              # dec 2014 meaning of pct_s2 changed
              pct_s2=100*sum(ntot>n1)/nsims, 
              nmean=mean(ntot), nrange=range(ntot), 
              nperc=quantile(ntot, p=npct))
  
  # table object summarizing the discrete distri of ntot
  # only if usePE=FALSE or if usePE=TRUE then a futility range must be used
  # or return it always?
  if (usePE==FALSE | (usePE==TRUE & fClower>0 & is.finite(fCupper))){
    res$ntable <- table(ntot)
  }
  
  if (details){
    cat("Total time (secs):\n")
    print(round((proc.time()-ptm),1))
    cat("\n")
  }
  
  class(res) <- c("pwrtsd", "list")
  return(res)
  
} #end function
