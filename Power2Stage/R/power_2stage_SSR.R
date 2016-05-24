# helper function to get CV from 'required sample size'
reqN2CV <- function(alpha=0.05, targetpower=0.8, theta0=1., theta1=0.8, n)
{
  z1 <- qnorm(1-alpha)
  if (log(theta0)!=0) z2 <- qnorm(targetpower) else z2 <- qnorm(1-(1-targetpower)/2) 

  s2e <- n*(log(theta0)-log(theta1))^2/2/(z1+z2)^2
  return(mse2CV(s2e))
}
  
# ----------------------------------------------------------------------------
# power (or alpha) of 2-stage studies with (blinded) sample size re-estimation
# see Golkowski et al.
#
# author D.L.
# ----------------------------------------------------------------------------
# require(PowerTOST)
# source("C:/Users/dlabes/workspace/Power2Stage/R/sampsiz2.R")
# source("C:/Users/dlabes/workspace/Power2Stage/R/sampsiz_n0.R")
# source("C:/Users/dlabes/workspace/Power2Stage/R/power.R")

power.2stage.ssr <- function(alpha=0.05, n1, GMR, CV, targetpower=0.8, 
                             pmethod=c("nct","exact", "shifted","ls"),
                             blind=FALSE, usePE=FALSE, min.n=0, max.n=Inf, 
                             theta0, theta1, theta2, npct=c(0.05, 0.5, 0.95), 
                             nsims, setseed=TRUE, details=FALSE)
{ # seems to give errors if alpha is a vector
  alpha <- alpha[1L]
  
  if (missing(CV)) stop("CV must be given!")
  if (CV<=0)       stop("CV must be >0!")
  
  if (missing(n1)) stop("Number of subjects in stage 1 must be given!")
  if (n1<=0)       stop("Number of subjects in stage 1 must be >0!")
  
  if (missing(GMR)) GMR <- 0.95
  
  if (missing(theta1) & missing(theta2))  theta1 <- 0.8
  if (!missing(theta1) & missing(theta2)) theta2 <- 1/theta1
  if (missing(theta1) & !missing(theta2)) theta1 <- 1/theta2
  
  if (GMR<=theta1 | GMR>=theta2) stop("GMR must be within acceptance range!")
  
  if (missing(theta0)) theta0 <- GMR
  
  if (n1>max.n) stop("max.n < n1 doestn't make sense!")
  if (min.n!=0 & min.n<n1) stop("min.n < n1 doestn't make sense!")
  
  if(missing(nsims)){
    nsims <- 1E5
    if(theta0<=theta1 | theta0>=theta2) nsims <- 1E6
  }
  
  # check if power calculation method is nct, exact, shifted or large sample
  pmethod <- match.arg(pmethod)
  
  if(blind & usePE) {
    usePE=FALSE
    warning("usePE=TRUE does not make sense for blinded SSR. Reset to usePE=FALSE.")
  }
  
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
  
# ----- interim ----------------------------------------------------------
# simulate blinded est. of variance
  Cfact <- bk/n1
  df    <- n1-2
  tval  <- qt(1-alpha[1], df)
  sdm   <- sqrt(mse*Cfact)
  # simulate point est. via normal distribution
  pes   <- rnorm(n=nsims, mean=mlog, sd=sdm)
  # simulate mse via chi-squared distribution
  mses  <- mse*rchisq(n=nsims, df=df)/df
   
  #s2os  <- mses + pes^2/2   # is this correct? example shows: it's not correct.
                             # s2os is smaller than mse! I would expect 
                             # the other way round!
  # SS of tmt formula is only valid if sequence balanced
  #        res. SS + SS of tmt 
  s2os <- (df*mses + n1*pes^2/2)/(n1-1)
  # unblind
  if (!blind) s2os <- mses

  if(details){
    # time for stage 1 sims
    cat(" - Time consumed (secs):\n")
    print(round((proc.time()-ptm),2))
    cat("Keep calm. ")
  }

  # ------ recalculate sample size -----------------------------------------
  ptms <- proc.time()
  # first calc. power to avoid unnecessary time consuming sample size estimation
  # this is not really part of the method (?) but a run-time goody
  # will give some boost for small CV's and/or high n1
  if(pmethod!="ls"){
    # use GMR or pe1 in sample size re-est.
    lpes <- lGMR
    if(usePE) lpes <- pes
    pwr <- .calc.power(alpha=alpha, ltheta1=ltheta1, ltheta2=ltheta2, 
                       diffm=lpes, sem=sqrt(bk*s2os/n1), df=df, method=pmethod)
    if(details){
      cat("Sample sizes (", sum(pwr<targetpower),
          " studies) will be re-estimated.\n", sep="")
      cat("May need some time.\n")
    }
  } else {
    # don't have power for large sample approx.
    # moreover ssr is so fast that we don't need the power step
    pwr  <- rep(0, times=nsims)
    if(details){
      cat("Sample sizes will be re-estimated.\n")
      cat("May need some time.\n")
    }
  }
  # total sample size
  ntot    <- rep(n1, times=nsims)
  mse_tmp <- s2os[pwr<targetpower]
  pes_tmp <- pes[pwr<targetpower]
  # use GMR or pe1 in sample size re-est.
  lpes <- lGMR
  if(usePE) lpes <- pes_tmp
  if(pmethod=="ls"){
    # large sample approx. via normal distribution
    nt <- .sampleN00(alpha=alpha, targetpower=targetpower, se=sqrt(mse_tmp), 
                     diffm=lpes, ltheta1=ltheta1, ltheta2=ltheta2, bk=2, 
                     steps=2, diffmthreshold=0.0)
  } else {
    # exact or approx. via t-distri
    nt <- .sampleN2(alpha=alpha, targetpower=targetpower, ltheta0=lpes,
                    mse=mse_tmp, ltheta1=ltheta1, ltheta2=ltheta2, 
                    method=pmethod, bk=2)
  }
  #browser()
  # maybe we have enough power in all cases and thus no re-estimated sample size
  if(length(nt)>0) ntot[pwr<targetpower] <- nt
  # take care of memory
  rm(nt, pwr, mse_tmp, lpes)
  # sample size returns Inf if pe outside acceptance range, then stay with n1 
  # but this should not occure here since pe1 is not used
  ntot <- ifelse(is.finite(ntot), ntot, n1)
  # use max.n if nt > max.n
  ntot <- ifelse(ntot>max.n, max.n, ntot)
  # use min.n if nt < min.n
  ntot <- ifelse(ntot<min.n, min.n, ntot)
  
  n2   <- ntot-n1
  # do not fall below n1
  n2 <- ifelse(n2<0, 0, n2)

  if(details){
    cat("Time consumed (secs):\n")
    print(round((proc.time()-ptms),1))
  }

  # ---------- stage 2 evaluation --------------------------------------
  m1    <- pes
  SS1   <- (n1-2)*mses
  # to avoid warnings for n2=0 in rnorm() and rchisq()
  ow    <- options("warn")
  options(warn=-1)
  m2    <- ifelse(n2>0, rnorm(n=nsims, mean=mlog, sd=sqrt(mse*bk/n2)), 0)
  # ??? (n2-2) cancels out! 
  SS2   <- ifelse(n2>2, (n2-2)*mse*rchisq(n=nsims, df=n2-2)/(n2-2), 0)
  # reset options
  options(ow)

  ntot <- n1+n2

  # do we need this 'stage' contribution?
  # can't discover such in Golkowski et al., but ... 
  withStage <- TRUE
  if(withStage){
    SSmean <- ifelse(n2>0, (m1-m2)^2/(2/n1+2/n2), 0)
    df2    <- ifelse(n2>0, ntot-3, n1-2)
  } else {
    SSmean <- 0
    df2    <- ntot-2
  }
  pe2    <- ifelse(n2>0, (n1*m1+n2*m2)/ntot, pes)
  mse2   <- ifelse(n2>0, (SS1+SSmean+SS2)/df2, mses)
  # take care of memory
  rm(m1, m2, SS1, SS2, SSmean)
  # calculate CI and decide BE 
  hw    <- qt(1-alpha, df2)*sqrt(mse2*bk/ntot)
  lower <- pe2 - hw
  upper <- pe2 + hw
  BE    <- (lower>=ltheta1) & (upper<=ltheta2)
  # take care of memory
  rm(lower, upper, hw)

  # the return list
  res <- list(method="SSR",
              alpha=alpha, CV=CV, n1=n1, GMR=GMR, targetpower=targetpower, 
              pmethod=pmethod, theta0=theta0, theta1=theta1, theta2=theta2, 
              usePE=usePE, max.n=max.n, min.n=min.n, blind=blind, nsims=nsims,
              # results 
              pBE=sum(BE)/nsims, 
              #pBE_s1 ?
              pct_s2=100*length(ntot[ntot>n1])/nsims, 
              nmean=mean(ntot), nrange=range(ntot), nperc=quantile(ntot, p=npct))
  
  # return a table object as summary of ntot distribution
  # only if usePE==FALSE ? or always?
  if(!usePE) res$ntable <- table(ntot)
  
  if (details){
    cat("Total time consumed (secs):\n")
    print(round((proc.time()-ptm),1))
    cat("\n")
  }
  
  # output is now via S3 print method

  class(res) <- c("pwrtsd", "list")
  return(res)
  
} #end function
