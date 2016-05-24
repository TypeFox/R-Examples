# --------------------------------------------------------------------------
# power (or alpha) of 2-stage studies with a 2-group parallel desig according 
# to Potvin et. al. # methods "B" and "C", modified to include a futility 
# criterion Nmax and modified to use PE of stage 1 in sample size estimation
#
# variant of power.2stage.p() wich calculates power always via pooled t-test
# formulas according to the fuglsang 2014 paper
#
# author D.L.
# --------------------------------------------------------------------------
# require(PowerTOST)

power.2stage.pAF <- function(method=c("B","C"), alpha0=0.05, alpha=c(0.0294,0.0294), 
                             n1, GMR, CV, targetpower=0.8, 
                             pmethod=c("shifted", "nct", "exact"), 
                             usePE=FALSE, Nmax=Inf, test=c("welch", "t-test", "anova"),
                             theta0, theta1, theta2, npct=c(0.05, 0.5, 0.95),  
                             nsims, setseed=TRUE, details=FALSE)
{
  if (missing(CV)) stop("CV(s) must be given!")
  if (any(CV<=0))  stop("CV(s) must be >0!")
  # equal CV's
  if (length(CV)==1) {
    CVT <- CVR <-CV
  } else {
    CV <- CV[1:2] # truncate if longer then 2
    CVT <- CV[1]; CVR <- CV[2]
  }
  varT <- CV2mse(CVT)
  varR <- CV2mse(CVR)
  
  if (missing(n1)) stop("Number of subjects in stage 1 must be given!")
  if (length(n1)>1) {
    warning("n1 has to be a scalar. Sum of vector will be used.")
    n1 <- sum(n1)
  }
  if (n1<=0)    stop("Number of subjects in stage 1 must be >0!")
  if (n1%%2!=0) warning("Number of subjects in stage 1 should be even.\n",
                        "  Will be truncated to even.", immediate. = TRUE)
  
  if (missing(GMR)) GMR <- 0.95
  
  if (missing(theta1) & missing(theta2))  theta1 <- 0.8
  if (!missing(theta1) & missing(theta2)) theta2 <- 1/theta1
  if (missing(theta1) & !missing(theta2)) theta1 <- 1/theta2
  
  if (GMR<=theta1 | GMR>=theta2) stop("GMR must be within acceptance range!")
  
  if (missing(theta0)) theta0 <- GMR
  
  if (n1>Nmax) stop("n1>Nmax doestn't make sense!")

  if(missing(nsims)){
    nsims <- 1E5
    if(theta0<=theta1 | theta0>=theta2) nsims <- 1E6
  }
  
  # check if Potvin B or C
  method  <- match.arg(method)
  # check test
  test  <- match.arg(test)
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
  lGMR    <- log(GMR)
  mlog    <- log(theta0)
  bk      <- 4   # 2-group parallel design constant
  # reserve memory for the BE result
  BE      <- rep.int(NA, times=nsims)
  
# ----- stage 1 ----------------------------------------------------------
  
  ns1T  <- ns1R <- trunc(n1/2) # if not even
  nT    <- ns1T; nR <- ns1R
  n1    <- ns1T + ns1R
  Cfact <- bk/n1
  df    <- ns1T + ns1R -2
  tval  <- qt(1-alpha[1], df) # ANOVA and t-test
  # simulate means via normal distributions
  m1T  <- rnorm(n=nsims, mean=mlog, sd=sqrt(varT/ns1T))
  m1R  <- rnorm(n=nsims, mean=0, sd=sqrt(varR/ns1R))
  # point est.
  pes   <- m1T-m1R
  # simulate variances via chi-squared distribution
  varsT  <- varT*rchisq(n=nsims, df=ns1T-1)/(ns1T-1)
  varsR  <- varR*rchisq(n=nsims, df=ns1R-1)/(ns1R-1)
  Vpooled <- ((ns1T-1)*varsT + (ns1R-1)*varsR)/(n1-2)

  if(method=="C"){
    # if method=C then calculate power for alpha0=0.05 and plan GMR
    # Here we use A.Fuglsangs settings, namely power monitoring steps
    # and sample size via pooled t-test
    pwr <- .calc.power(alpha=alpha0, ltheta1=ltheta1, ltheta2=ltheta2, 
                       diffm=lGMR, sem=sqrt(bk*Vpooled/n1), df=df,  
                       method=pmethod)
    # calculate CIs at alpha0 for all
    if (test=="t-test" || test=="anova"){
      tval0 <- qt(1-alpha0, df)
      hw    <- tval0*sqrt(bk*Vpooled/n1)
    } else {
      # Welch's t-test
      se  <- sqrt(varsT/nT + varsR/nR)
      dfs <- (varsT/nT + varsR/nR)^2/(varsT^2/nT^2/(nT-1)+varsR^2/nR^2/(nR-1))
      # dfs <- trunc(dfs)
      hw  <- qt(1-alpha0,dfs)*se
    }
    lower <- pes - hw
    upper <- pes + hw
    # fail or pass
    BE    <- lower>=ltheta1 & upper<=ltheta2
    # if power>0.8 then calculate CI for alpha=0.05
    # i.e. if power<0.8 then 
    BE[pwr<targetpower] <- NA # not yet decided
  }
  # method "B" or power<=0.8 in method "C"
  # evaluate BE with alpha[1]
  Vpooled_tmp <- Vpooled[is.na(BE)]
  pes_tmp <- pes[is.na(BE)]
  varsT_tmp   <- varsT[is.na(BE)] 
  varsR_tmp   <- varsR[is.na(BE)] 
  BE1 <- rep.int(NA, times=length(Vpooled_tmp))
  # calculate CI for alpha=alpha1
  if (test=="t-test" || test=="anova"){
    hw    <- tval*sqrt(bk*Vpooled_tmp/n1)
  } else {
    # Welch's t-test
    se  <- sqrt(varsT_tmp/nT + varsR_tmp/nR)
    dfs <- (varsT_tmp/nT + varsR_tmp/nR)^2/(varsT_tmp^2/nT^2/(nT-1) + 
                                    varsR_tmp^2/nR^2/(nR-1))
    # dfs <- trunc(dfs)
    hw  <- qt(1-alpha[1],dfs)*se
    rm(se, dfs)
  }
  rm(varsT_tmp, varsR_tmp)
  lower <- pes_tmp - hw
  upper <- pes_tmp + hw
  BE1   <- lower>=ltheta1 & upper<=ltheta2
  if (method=="C"){
    #if BE met -> PASS and stop
    #if not BE -> goto sample size estimation i.e flag BE1 as NA
    BE1[!BE1] <- NA
  } else { 
    # method B
    # evaluate power at alpha[1] and planGMR
    pwr <- .calc.power(alpha=alpha[1], ltheta1=ltheta1, ltheta2=ltheta2, 
                       diffm=lGMR, sem=sqrt(bk*Vpooled_tmp/n1), df=df, 
                       method=pmethod)
    # if BE met then decide BE regardless of power
    # if not BE and power<0.8 then goto stage 2
    BE1[ !BE1 & pwr<targetpower] <- NA 
  }
  # combine 'stage 0' from method C and stage 1
  BE[is.na(BE)] <- BE1
  # take care of memory
  # done with them
  rm(BE1, hw, lower, upper, pes_tmp, Vpooled_tmp)
  
  # time for stage 1
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
      cat("Keep calm. Sample sizes for stage 2 (", length(pes_tmp),
          " studies)\n", sep="")
      cat("will be estimated. May need some time.\n")
    }
    # preliminary setting stage=2 for those not yet decided BE
    # may be altered for those with nts>Nmax or nts=Inf 
    # from sample size est. if pe outside acceptance range
    # see below
    stage[is.na(BE)] <- 2
    Vpooled_tmp <- Vpooled[is.na(BE)]
    m1T <- m1T[is.na(BE)]
    m1R <- m1R[is.na(BE)]
    varsT <- varsT[is.na(BE)]
    varsR <- varsR[is.na(BE)]
    BE2      <- rep.int(NA, times=length(Vpooled_tmp))
    s2       <- rep.int( 2, times=length(Vpooled_tmp))
    #------ sample size for stage 2 ---------------------------------------
    ptms <- proc.time()
    if (usePE){
      # use mse1 & pe1 like in the paper of Karalis/Macheras
      # sample size function returns Inf if pe1 is outside acceptance range
#       nts <- mapply(FUN=.sampleN, mse=Vpooled_tmp, ltheta0=pes_tmp, 
#                     MoreArgs=list(alpha=alpha[2], targetpower=targetpower, 
#                                   ltheta1=ltheta1, ltheta2=ltheta2,
#                                   method=pmethod, bk=4))
      nts <- .sampleN2(alpha=alpha[2], targetpower=targetpower, ltheta0=pes_tmp,
                       mse=Vpooled_tmp, ltheta1=ltheta1, ltheta2=ltheta2, 
                       method=pmethod, bk=4)
    } else {
      # use mse1 & plan GMR to calculate sample size (original Potvin)
#       nts <- mapply(FUN=.sampleN, mse=Vpooled_tmp, 
#                     MoreArgs=list(alpha=alpha[2], targetpower=targetpower, 
#                                   ltheta0=lGMR, ltheta1=ltheta1, ltheta2=ltheta2,
#                                   method=pmethod, bk=4))
      nts <- .sampleN2(alpha=alpha[2], targetpower=targetpower, ltheta0=lGMR,
                       mse=Vpooled_tmp, ltheta1=ltheta1, ltheta2=ltheta2, 
                       method=pmethod, bk=4)
    }
    
    n2  <- ifelse(nts>n1, nts - n1, 0)
    
    if(details){
      cat("Time consumed (secs):\n")
      print(round((proc.time()-ptms),1))
    }
    # futility rule: if nts > Nmax -> stay with stage 1 result not BE
    # ntotal = n1 reasonable?
    if (is.finite(Nmax) | any(!is.finite(nts))){
      # sample size may return Inf if PE is used in ss estimation
      # in that case we stay with stage 1
      BE2[!is.finite(n2) | (n1+n2)>Nmax] <- FALSE
      # and we are counting these for stage 1
      s2[BE2==FALSE]  <- 1
      # debug print
      # cat(sum(!BE2, na.rm=T)," cases with nts>Nmax or nts=Inf\n")
      # save 
      stage[is.na(BE)] <- s2
      # save the FALSE and NA in BE
      BE[is.na(BE)]    <- BE2
      # filter out those were BE was yet not decided
      pes_tmp  <- pes_tmp[is.na(BE2)]
      m1T <- m1T[is.na(BE2)]
      m1R <- m1R[is.na(BE2)]
      Vpooled_tmp <- Vpooled_tmp[is.na(BE2)]
      varsT <- varsT[is.na(BE2)]
      varsR <- varsR[is.na(BE2)]
      n2    <- n2[is.na(BE2)]
    } # end of futility Nmax
    # ----- simulate stage 2 data ------
    nsim2 <- length(pes_tmp)
    ns2T  <- ns2R <- n2/2
    # to avoid warnings for ns2X==0 in rnorm() and ns2X<=1 in rchisq()
    ow    <- options("warn")
    options(warn=-1)
    ms2T  <- ifelse(ns2T>0, rnorm(n=nsim2, mean=mlog, sd=sqrt(varT/ns2T)), 0)
    ms2R  <- ifelse(ns2R>0, rnorm(n=nsim2, mean=0, sd=sqrt(varR/ns2R)), 0)
    # means T/R
    nT <- ns1T+ns2T
    nR <- ns1R+ns2R
    mT <- (ns1T*m1T+ns2T*ms2T)/nT
    mR <- (ns1R*m1R+ns2R*ms2R)/nR
    # point est.
    pes <- mT-mR
    # means for s1, s2 (stages)
    m1  <- (ns1T*m1T + ns1R*m1R)/(ns1T+ns1R)
    m2  <- ifelse((ns2T+ns2R)>0, (ns2T*ms2T + ns2R*ms2R)/(ns2T+ns2R),0)
    # total mean
    mt  <- (nT*mT+nR*mR)/(nT+nR)
    # simulate variances via chi-squared distribution
    # attention! in case of ns2X==1 rchisq gives NaN!
    # TODO: work out the correct way for ns2X==1 
    vars2T  <- ifelse(ns2T>1, varT*rchisq(n=nsim2, df=ns2T-1)/(ns2T-1), 0)
    vars2R  <- ifelse(ns2R>1, varR*rchisq(n=nsim2, df=ns2R-1)/(ns2R-1), 0)
    # reset options
    options(ow) 
    
    # vars T/R over stage 1, stage 2
    # s2y = sum of y squared
    # if ns2X==0 then we get here originally NA
    #sy2T <- (ns1T-1)*varsT+m1T^2/ns1T + (ns2T-1)*vars2T+ms2T^2/ns2T
    sy2T <- (ns1T-1)*varsT+m1T^2/ns1T
    sy2T <- ifelse(ns2T>0, sy2T+(ns2T-1)*vars2T+ms2T^2/ns2T, sy2T)
    sy2R <- (ns1R-1)*varsR+m1R^2/ns1R
    sy2R <- ifelse(ns2R>0, sy2R+(ns2R-1)*vars2R+ms2R^2/ns2R, sy2R)
    varsT <- (sy2T - mT^2/nT)/(nT-1)
    varsR <- (sy2R - mR^2/nR)/(nR-1)
    rm(sy2T, sy2R)
    # calculate CI for stage 2 with alpha[2]
    if (test=="t-test" || test=="anova"){
      Vpooled <- ((nT-1)*varsT + (nR-1)*varsR)/(nT+nR-2)
      dfs <- (nT+nR-2)
      if (test=="anova"){
        # no subjects in stage 2 ((ns2R+ns2T)==0) may occure in case of 
        # high n1 and/or haybittle-peto alpha's
        Vpooled <- ifelse(ns2R+ns2T>0,
                         (dfs*Vpooled - n1*(m1-mt)^2 - n2*(m2-mt)^2)/(dfs-1),
                          Vpooled)
        dfs     <- ifelse(ns2R+ns2T>0, dfs-1, dfs)
      }
      hw    <- qt(1-alpha[2],dfs)*sqrt(bk*Vpooled/(nT+nR))
      rm(Vpooled, dfs)
    } else {
      # Welch's t-test
      se  <- sqrt(varsT/nT + varsR/nR)
      dfs <- (varsT/nT + varsR/nR)^2/(varsT^2/nT^2/(nT-1)+varsR^2/nR^2/(nR-1))
      #dfs <- trunc(dfs)
      hw  <- qt(1-alpha[2],dfs)*se
      rm(se, dfs)
    }
    
    lower <- pes - hw
    upper <- pes + hw
    BE2   <- lower>=ltheta1 & upper<=ltheta2
    # combine stage 1 & stage 2
    ntot[is.na(BE)]  <- n1+n2
    BE[is.na(BE)]    <- BE2
    # done with them
    rm(BE2, nts, lower, upper, hw, ns2T, ns2R, ms2T, ms2R)
  } # end stage 2 calculations
  # take care of memory
  rm(pes_tmp, pes)
  # the return list
  res <- list(design="2 parallel groups", method=method, 
              alpha0=ifelse(method=="C",alpha0,NA), alpha=alpha, CV=CV, n1=n1, 
              GMR=GMR, test=test, targetpower=targetpower, pmethod=pmethod, 
              theta0=exp(mlog), theta1=theta1, theta2=theta2, 
              usePE=usePE, Nmax=Nmax, nsims=nsims,
              # results 
              pBE=sum(BE)/nsims, pBE_s1=sum(BE[stage==1])/nsims,
              # Dec 2014: meaning of pct_s2 changed
              pct_s2=100*sum(ntot>n1)/nsims, 
              nmean=mean(ntot), nrange=range(ntot), nperc=quantile(ntot, p=npct))
  
  # table object summarizing the discrete distri of ntot
  # only if usePE=FALSE or if usePE=TRUE then Nmax must be finite
  # or return it always?
  if (usePE==FALSE | (usePE==TRUE & is.finite(Nmax))){
    res$ntable <- table(ntot)
  }
  if (details){
    cat("Total time consumed (secs):\n")
    print(round((proc.time()-ptm),1))
    cat("\n")
  }
    
  class(res) <- c("pwrtsd", "list")
  return(res)
  
} #end function
