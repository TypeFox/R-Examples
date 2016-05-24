# --------------------------------------------------------------------------
# power (or alpha) of 2-stage studies according to Potvin et. al. 
# methods "B" and "C", 
# - modified to include a futility criterion Nmax
# - modified to use PE of stage 1 in sample size estimation if PE1 is
#   outside GMR and 1/GMR
#
# Author D.L.
# --------------------------------------------------------------------------
# require(PowerTOST)
# source("./R/sampsiz2.R")
# source("./R/sampsiz_n0.R")
# source("./R/power.R")

power.2stage.DL <- function(method=c("B","C"), alpha0=0.05, alpha=c(0.0294,0.0294),
                         n1, GMR, CV, targetpower=0.8, 
                         pmethod=c("nct","exact", "shifted"),
                         Nmax=Inf, min.n2=0, theta0, theta1, theta2,  
                         npct=c(0.05, 0.5, 0.95), nsims, setseed=TRUE, 
                         print=TRUE, details=FALSE)
{
  if (missing(CV)) stop("CV must be given.")
  if (CV<=0)       stop("CV must be >0.")
  
  if (missing(n1)) stop("Number of subjects in stage 1 must be given.")
  if (n1<=0)       stop("Number of subjects in stage 1 must be >0.")
  
  if (missing(GMR)) GMR <- 0.95
  
  if (missing(theta1) & missing(theta2))  theta1 <- 0.8
  if (!missing(theta1) & missing(theta2)) theta2 <- 1/theta1
  if (missing(theta1) & !missing(theta2)) theta1 <- 1/theta2
  
  if (GMR<=theta1 | GMR>=theta2) stop("GMR must be within acceptance range.")
  
  if (missing(theta0)) theta0 <- GMR
  
  if (n1>Nmax) stop("n1>Nmax doesn't make sense!")
  
  if(missing(nsims)){
    nsims <- 1E5
    if(theta0<=theta1 | theta0>=theta2) nsims <- 1E6
  }
  
  if(min.n2!=0 & min.n2<2) stop("min.n2 has to be at least +2.")
  # make even (round up)
  if( min.n2%%2 != 0) {
    min.n2 <- min.n2 + min.n2%%2
    message("min.n2 rounded to even", min.n2)
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
  sdm   <- sqrt(mse*Cfact)
  # simulate point est. via normal distribution
  pes   <- rnorm(n=nsims, mean=mlog, sd=sdm)
  # simulate mse via chi-squared distribution
  mses  <- mse*rchisq(n=nsims, df=df)/df
  
  if(method=="C"){
    # if method=C then calculate power for alpha0=0.05 and plan GMR
    pwr <- .calc.power(alpha=alpha0, ltheta1=ltheta1, ltheta2=ltheta2, 
                       diffm=lGMR, sem=sqrt(bk*mses/n1), df=df, method=pmethod)
    
    tval0 <- qt(1-alpha0, df)
    hw    <- tval0*sqrt(Cfact*mses)
    lower <- pes - hw
    upper <- pes + hw
    # fail or pass
    BE    <- lower>=ltheta1 & upper<=ltheta2
    # if power>0.8 then calculate CI for alpha=0.05
    # i.e. if power<0.8 then 
    BE[pwr<targetpower] <- NA # not yet decided
  }
  # method "B" or power<=0.8 in method "C"
  # calculate power for alpha=alpha[1]
  mses_tmp <- mses[is.na(BE)]
  pes_tmp  <- pes[is.na(BE)]
  BE1 <- rep.int(NA, times=length(mses_tmp))
  # calculate CI for alpha=alpha1
  hw    <- tval*sqrt(Cfact*mses_tmp)
  lower <- pes_tmp - hw
  upper <- pes_tmp + hw
  BE1   <- lower>=ltheta1 & upper<=ltheta2
  if (method=="C"){
    #if BE met -> PASS stop
    #if not BE -> goto sample size estimation i.e flag BE1 as NA
    BE1[!BE1] <- NA
  } else { 
    # method B
    # evaluate power at alpha[1]
    pwr <- .calc.power(alpha=alpha[1], ltheta1=ltheta1, ltheta2=ltheta2, 
                       diffm=lGMR, sem=sqrt(bk*mses_tmp/n1),df=df,  
                       method=pmethod)
    # if BE met then decide BE regardless of power
    # if not BE and power<0.8 then goto stage 2
    BE1[ !BE1 & pwr<targetpower] <- NA 
  }
  # combine 'stage 0' from method C and stage 1
  BE[is.na(BE)] <- BE1
  # take care of memory
  # done with them
  rm(BE1, hw, lower, upper)
  
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
    # may be altered for those with nt>Nmax or nt=Inf 
    # from sample size est. if pe outside acceptance range
    # see below
    stage[is.na(BE)] <- 2
    mses_tmp <- mses[is.na(BE)]
    BE2      <- rep.int(NA, times=length(mses_tmp))
    s2       <- rep.int(2, times=length(mses_tmp))
    #------ sample size for stage 2 ---------------------------------------
    ptms <- proc.time()
    # use mse1 & pe1 or mse1 & GMR
    # sample size function returns Inf if pe1 is outside acceptance range
    #browser()
    # now the rule: PE1 in range GMR ... 1/GMR then use GMR, else use PE1
    # pes_tmp is needed later, thus create another variable
    lGMR2 <- ifelse(lGMR<0, lGMR, -lGMR)
    pes_tmp2 <- ifelse(pes_tmp>=lGMR2 & pes_tmp<=-lGMR2, lGMR, pes_tmp )
    nt <- .sampleN2(alpha=alpha[2], targetpower=targetpower, ltheta0=pes_tmp2,
                    mse=mses_tmp, ltheta1=ltheta1, ltheta2=ltheta2, method=pmethod)
    rm(pes_tmp2)
    
    n2 <- ifelse(nt>n1, nt - n1, 0)
    # assure a min.n2
    n2 <- ifelse(n2<min.n2, min.n2, n2)
    # some upper bound due to numerics?
    #n2  <- ifelse(n2>100000, 100000, n2) # may not necessary
    
    if(details){
      cat("Time consumed (secs):\n")
      print(round((proc.time()-ptms),1))
    }
    # futility rule: if nt > Nmax -> stay with stage 1 result: not BE
    if (is.finite(Nmax) | any(!is.finite(nt))){
      # sample size may return Inf if PE is used in ss estimation
      # in that case we stay with stage 1
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
    # ??? (n2-2) cancels out! 
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
    #browser()
    # combine those with stage 1 only & those having stage 2
    ntot[is.na(BE)]  <- nt
    BE[is.na(BE)]    <- BE2
    # done with them
    rm(BE2, nt, lower, upper, hw)
  } # end stage 2 calculations
  # take care of memory
  rm(pes_tmp, mses_tmp)
  # the return list
  res <- list(design="2x2 crossover",
              method=method, modified="DL",
              alpha0=ifelse(method=="C",alpha0,NA), alpha=alpha, 
              CV=CV, n1=n1, GMR=GMR, targetpower=targetpower, pmethod=pmethod, 
              theta0=exp(mlog), theta1=theta1, theta2=theta2, usePE=NULL, 
              Nmax=Nmax, min.n2=min.n2, nsims=nsims,
              # results 
              pBE=sum(BE)/nsims,
              # Dec 2014 changed to
              pBE_s1=sum(BE[ntot==n1])/nsims,
              # was
              #pBE_s1=sum(BE[stage==2])/nsims,
              # Dec 2014 changed the meaning of pct_s2, was
              #pct_s2=100*length(BE[stage==2])/nsims,
              # whereby in case of unsymmetric alpha's stage 2 was also assigned
              # if n2=0 but not BE using alpha[1] to fascilate BE test with alpha[2]
              # now it is 
              pct_s2=100*sum(ntot>n1)/nsims, 
              # which simply means all those with n2>0
              nmean=mean(ntot), nrange=range(ntot), nperc=quantile(ntot, p=npct)
              )

  # table object summarizing the discrete distri of ntot
  # here usePE=TRUE, then Nmax must be finite?
  if (is.finite(Nmax))   res$ntable <- table(ntot)
  #res$ntable <- table(ntot)
  
  if (details){
    cat("Total time consumed (secs):\n")
    print(round((proc.time()-ptm),1))
    cat("\n")
  }
  
  class(res) <- c("pwrtsd", "list")
  return(res)
  
} #end function
