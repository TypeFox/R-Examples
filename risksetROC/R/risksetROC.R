## REFERENCE: HAEGERTY & ZHENG, 2005, Biometrics 61, 92-105
## THIS FILE CONTAINS ALL THE NECESARRY R-CODES FOR CREATING
## RISKSETROC (INCIDENT/DYNAMIC) AS DISCUSSED IN THE PAPER
## ABOVE.
## AUTHOR/COMPILER: P. SAHA
## R-CODES ORIGINALLY WRITTEN BY: P. HEAGERTY (FUNCTIONS 1-7)
## DATE: AUG 17, 2006
## ADDITIONAL R-CODES CREATED BY: P. SAHA (FUNCTIONS 8-9)
## DATE: OCT 13, 2006

## LIST OF FUNCTIONS AND THEIR PURPOSE:

## This file contains the following functions:
## 1. CoxWeights(eta, Stime, status, predict.time)
## 2. IntegrateAUC(AUC, utimes, St, tmax, weight="rescale")
## 3. llCoxReg3(Stime, status, x, span=0.40, p=1, window="asymmetric", shrink=F)
## 4. riskset(dat)
## 5. weightedKM(time, Sstatus, wt=NULL, time2=NULL)
## 6. KM.plot(times, survival, max.T=NULL, lty=NULL, all=T)
## 7. SchoenSmooth(fit, Stime, status, span=0.40, order=0, shrink=F)
## Description of the functions preecede the function definitions
##
## Date: Oct 11, 2006
## We include two functions risksetROC() and risksetAUC() which uses the
## other functions  implicitely rather than calling the functions on
## their own.

## 1. CoxWeights(marker, time, status, predict.time, time2=NULL):
## Returns estimates of TP and FP based on a Cox PH model. TP is
## estimated as Equation (1) and FP is estimated as Equation (2)
## of the paper. Here eta is the estimated linear predictor as
## obtained from
## a Cox model. The other arguments are time denoting the survival
## time or entry time, time2 denoting the exit time if any and
## status denoting the status of the patients (1: dead,
## 0:censored). predict.time denotes the time point of interest.

## #############################
## THIS FUNCTION IS WORKING A OK
## #############################

CoxWeights <- function(marker, Stime, status, predict.time, entry=NULL)
  ## DATE: OCT 18, 2006 -- argument eta changed to marker
  ## DATE: NOV 20, 2006 -- argument Stime changed to time
  ##                       argument time2 introduced
  ## DATE: Feb 12, 2007 -- argument time changed to Stime
  ##                       argument time2 changed to entry
{
  eta <- marker
  Target <- predict.time
  if(length(entry) == 0)
    {
      entry=rep(0,NROW(Stime))
    }

  at.risk <- ((Stime>=Target)&(entry<=Target))
  
  the.eta <- eta[ at.risk ] ## markers of the subjects in the riskset
  n <- length(the.eta) ## how many are in the riskset 
  the.dead <- (Stime==Target)&(status==1) ## how many died at predict.time
  the.dead <- the.dead[ at.risk ]
  n.dead <- sum( the.dead )
  p0 <- rep( 1/(n-n.dead), n )
  if( n.dead>0 ) p0[ the.dead ] <- 0
  p1 <- exp( the.eta )/sum( exp(the.eta) )
  ooo <- order(the.eta)
  eta.out <- the.eta[ ooo ]
  TP <- c( 1, 1 - cumsum( p1[ ooo ] ), 0 )
  FP <- c( 1, 1 - cumsum( p0[ ooo ] ), 0 )

  dFP <- abs( FP[-1] - FP[-length(FP)] )
  aTP <- 0.5*( TP[-1] + TP[-length(TP)] )
  area <- sum( dFP * aTP )

  out <- list( marker = eta.out, TP=TP, FP=FP, AUC=area )
  return(out)
}



## 2. IntegrateAUC( AUC, utimes, St, tmax, weight="rescale" )
## This function integrates AUC using w(t) = 2*f(t)*S(t).
## The arguments are AUC at ordered unique failure times utimes.
## Survival estimates St at those utimes. tmax denotes the max
## time length we want to consider. Two options for the argument
## weight is given. These options are "rescale" and  "conditional". 
## If weight = "rescale" then the weights are rescaled so that the
## sum of the weights is unity. If weight = "conditional" then the
## weights are calculated conditioned on the fact that max time is
## some finite (given) time.
##*** change the order of AUC and utimes?

## #############################
## THIS FUNCTION IS WORKING A OK
## #############################

IntegrateAUC <- function( AUC, utimes, St, tmax, weight="rescale" )
{
  ##
  ## Assume ordered utimes
  ##
  wChoice <- match( weight, c("rescale","conditional") )

  if( is.na(wChoice) )
    {
    cat("error in weight choice")
    stop(0)
  }
  ft <- rep( NA, length(St) )
  ft[1] <- 1.0 - St[1]
  ## this is original line:
  ## for( j in 1:length(St) ) ft[j] <- St[j-1] - St[j]
  ## this is corrected line: (?)
  for( j in 2:length(St) ) ft[j] <- St[j-1] - St[j]
  ##
  mIndex <- length( utimes[ utimes <= tmax ] )
  www <- 2*ft*St
  iAUC <-  0.0
  wTotal <- sum( www[1:mIndex] )
  Smax <- St[ min(mIndex+1,length(St)) ]
  ##   for( j in 1:mIndex ){
  ##     if( wChoice==1 ){
  ##       wj <- 2*ft[j]*St[j]/wTotal
  ##     }
  ##     if( wChoice==2 ){
  ##       wj <- 2*( ft[j]/(1-Smax) )*(St[j]-Smax)/(1-Smax)
  ##     }
  ##     iAUC <- iAUC + wj*AUC[j]
  ##   }
  ##   iAUC
  ##   
  if(wChoice == 1)      w=2*ft*St/wTotal
  else  w = 2*ft*(St-Smax)/((1-Smax)^2)
  iAUC=sum(w[1:mIndex]*AUC[1:mIndex])
  return(iAUC)
}

## *** NOTE:*** Example VA data
## Time taken with for loop: 0.002 0.000 0.001 0.000 0.000
## Time taken without for loop: 0.001 0.000 0.001 0.000 0.000
## Also the values computed with for loop and w/o it are the same,
## so taking out the for loop.

## 3. llCoxReg(time, time2=NULL, status, x, span=0.40, p=1,
##   window="asymmetric") :
## Calculates the parameter estimate \beta(t) of non-proportional
## hazard model using local-linear Cox regression. Also includes a
## function named riskset() to construct risk set at each  unique
## failure time.  Uses the riskset() function given afterwards.
## For p=1, the first column of beta gives the estimate of time-varying
## parameter and the second column gives the estimate of derivative
## of the time-varying  parameter. However, if the derivative of the
## time-varying covariate is of interest, then we suggest using p=2,
## where the first column of beta estimates the time-varying
## covariate andd the second estimates its derivative.
## *** Oct 6,2006: argument shrink=F/T is removed

## #############################
## THIS FUNCTION IS WORKING A OK
## #############################

llCoxReg <- function(Stime, entry=NULL, status, marker, span=0.40, p=1,
                     window="asymmetric") 
{
  ##
  ## time-dependent version of cox model (Cai and Sun, 2003)
  ##
  x <- marker
  shrink=FALSE

  
  TheTime=Stime
  TheStatus <- status
  survival.status <- TheStatus
  survival.time <- TheTime
  eta <- x
  utimes <- unique( TheTime[ TheStatus==1 ] )
  ## unique exit times
  utimes <- utimes[ order(utimes) ]
  ##
  ## find the bandwidth for each time
  ##
  grid.t <- utimes
  nt <- length(grid.t)
  ##
  if( window=="asymmetric" )
    {
      lambda <- matrix(0, nt, 2)
      nd <- round(span*nt/2)
      ## nd*2 is the length of the nhbd 
      for (i in 1:nt)
        {
          delta <-  grid.t-grid.t[i] 
          keep <- delta <= 0
          if( sum(keep) > nd )
            {
              lambda[i,1] <- sort( abs(delta[keep]) )[nd]
            }
          else
            {
              lambda[i,1] <- max( abs(delta[keep]) )
            }
          ## if there is more than nd in the lower nhbd of the unique
          ## failure time then order the abs(diffs) and take  nd-th time
          ## of them, o/w consider all and take the max of abs(diffs)
          
          keep <- delta >= 0
          if( sum(keep) > nd )
            {
              lambda[i,2] <- sort( abs(delta[keep]) )[nd]
            }
          else
            {
              lambda[i,2] <- max( abs(delta[keep]) )
            }

          ## if there is more than nd in the upper nhbd of the unique
          ## failure time then order the abs(diffs) and take  nd-th time
          ## of them, o/w consider all and take the max of abs(diffs)
          ## lambda consists of the end points of nhbd intervals

        }  
    }
  else
    {
      lambda <- rep(0,nt)
      nd <- round( span*nt )
      ##
      for (i in 1:nt)
        {
          lambda[i] <- sort( abs( grid.t-grid.t[i] ) )[nd]
        }
    }
  ##
  ##
  ##
  if(length(entry)==0)
    {
      entry=rep(0, NROW(time))
    }
  tp <-cbind(entry, Stime, status, eta )

  ooo <- order(survival.time)
  tp <- tp[ooo,]
  ## order the data according to increasing exit time
  ##
  ## create time-dependent covariate
  ##

  ## ##########
  ## cat(str(tp))
  ## ##########
  
##   if(length(time2)==0)
##     {
##       tpp <- riskset(tp, time2=FALSE)
## 
##       ## ##########
##       cat("If LOOP")
##       ## ##########      
##     }
##   else
##     {
##
##  tpp <- riskset(tp, time2=TRUE )
## 
##       ## ##########
##       cat("ELSE LOOP")
##       ## ##########      
## 
##    }


  tpp <- riskset(tp, entry=TRUE )
  
  ## ##########
  ## cat(str(tpp))
  ## ##########

  ## tpp <- riskset(tp, time2=TRUE )
  ## creates riskset from tp for each unique failure time and adds first
  ## two columns as the start and end of intervals. Also removes the
  ## Stime  and status column from tp and adds a newstatus col.
  ## details below in the function definition

  nc <- ncol(tpp)
  ## ncol for tpp, same as ncol for tp (ncol(tp) - 2 (remove Stime and
  ## status) + 3 (adds risk set interval endpoints ane new status) )
  ## so if tp had 3 cols (Stime, status, x), the tpp would have 4 cols:
  ## start, finish, newstatus and x
  nu <- length(utimes)

  ## p denotes the power, two choices for p: 1 or 2.
  ## p = 1 fits a (local) Cox model with covariates:
  ## x and x*(ul of time interval - unique failure time (say j) in that
  ## interval) (=x*(tpp[,2]-utime[j]))
  ## p = 2 fits a (local) quadratic Cox model with the covariates: x,
  ## x*(tpp[,2]-utime[j]),   x^2, x^2 * (tpp[,2]-utime[j])
  ## the local Cox model is fitted at each unique failure times.

  if( p==1 )
    {
      betat <- matrix(0,nu,2)
    }
  if( p==2 )
    {
      betat <- matrix(0,nu,4)
    }
  ##
  shrinkage <- rep( NA, nu )
  ##
  for (j in 1:nu) {
    ##
    tt <- grid.t[j]
    ##
    ## use Epanechnikov's optimal kernel
    ##
    if( window=="asymmetric" )
      {
        lLow <- max( lambda[j,1], 1e-6 )
        lHigh <- max( lambda[j,2], 1e-6 )
        u <- ( (tpp[,2]-tt)/lLow )*( tpp[,2] <= tt ) +
          ( (tpp[,2]-tt)/lHigh )*( tpp[,2] > tt ) 
        iu <- ifelse( abs(u)>1 , 0, 1 )
        wt <- iu*0.75*(1-u^2)
      }
    else
      {
        u <- (tpp[,2]-tt)/lambda[j]
        iu <- ifelse(abs(u)>1,0,1)
        wt <- iu*0.75*(1-u^2)
      }
    ## no 0 wt allowed in cox model
    wt<-ifelse(wt==0,0.0000001,wt)
    ##
    if( p==1 )
      {
        tpcox  <- list( start  = tpp[,1],
                       finish = tpp[,2],
                       status = tpp[,3],
                       tpx    = tpp[,4],
                       tpxt   = tpp[,4]*(tpp[,2]-tt),
                       wt=wt )
        ##
        tpfit <- coxph( Surv(start,finish,status) ~ tpx+tpxt, tpcox,
                       weights=wt)
        ddd <- 2*abs( tpfit$loglik[1] - tpfit$loglik[2] )
        ddd <- (ddd - length(tpfit$coef))/ddd
        ##
        betat[j,] <- tpfit$coef[1:2]
        if( shrink ) betat[j,] <- ddd * tpfit$coef[1:2]
        shrinkage[j] <- ddd
      }
    if( p==2 ){
      ## tpcox  <- list( start  = tpp[,1],
      ##                finish = tpp[,2],
      ##                status = tpp[,3],
      ##                tpx    = tpp[,4],
      ##                tpxt   = tpp[,4]*(tpp[,2]-tt),
      ##                tpx2   = tpp[,4]^2,
      ##                tpx2t   = (tpp[,4]^2)*(tpp[,2]-tt),
      ##                wt=wt )
      ##
      ## tpfit <- coxph( Surv(start,finish,status) ~ tpx + tpxt +
      ##                tpx2 + tpx2t,
      ##                tpcox,
      ##                weights=wt)
      ## DATE:OCT 12, 2006 -- The above definition of tpcox and
      ## consequently tpfit changed to the following 
      tpcox <- list(start  = tpp[,1],                 
                    finish = tpp[,2],                  
                    status = tpp[,3],                  
                    tpx    = tpp[,4],                  
                    tpxt   = tpp[,4]*(tpp[,2]-tt),     
                    tpxt2  = tpp[,4]*((tpp[,2]-tt)^2),
                    wt=wt)
      tpfit <- coxph( Surv(start,finish,status) ~ tpx + tpxt +
                     tpxt2,
                     tpcox,
                     weights=wt)
      
      ddd <- 2*abs( tpfit$loglik[1] - tpfit$loglik[2] )
      ddd <- (ddd - length(tpfit$coef))/ddd
      ##
      betat[j,] <- tpfit$coef[1:3]
      if( shrink ) betat[j,] <- ddd * tpfit$coef[1:3]
      shrinkage[j] <- ddd
    }        
  }
  ##
  ## out <- list( time=utimes, beta=betat, shrinkage=shrinkage )
  ## DATE: OCT 18, 2006 -- shrinkage REMOVED FROM RETURN LIST
  out <- list( time=utimes, beta=betat)
  return(out)
}




## 4.riskset(dat): creates risk set at each unique failure time.
## DATE: NOV 17, 2006
## AUTHOR: PARAMITA SAHA
## NEW riskset FUNCTION WHICH ALLOWS ENTRY TIME
## DATA SHOULD CONTAIN ENTRY TIME OR TIME (IN CASE ONLY RIGHT CENSORED
## DATA), EXIT TIME(NULL IN CASE OF RIGHT CENSORED DATA), SURVIVAL
## STATUS AND OTHER COVARIATES 
## first sort the data according to entry time

## #############################
## THIS FUNCTION IS WORKING A OK
## #############################

riskset<-function(dat, entry=FALSE)
  ## creates risk set at each unique failure time
{
  ## tp has at least three variables: survival.times, survival.status and
  ## marker x.
  ## note that the data need not be ordered wrt increasing
  ## survival.times
  ## order the survival times if exit==NULL and call riskset()
  time2=entry

  tChoice <- match(time2, c("TRUE","FALSE"))  

  if(tChoice==2)
    {
      nobs=NROW(dat)
      dat=cbind(rep(0,nobs), dat)
    }

  ## cat(str(dat))

  time=dat[,1]
  time2=dat[,2]
  status=dat[,3]
  time2[status==1]=ifelse(
         time2[status==1]-time[status==1]==0,
         time2[status==1]+0.0001,
         time2[status==1])
  ## if a subject has entered and exited at the same time
  ## with an event, then add 0.0001 to his exit time      
  dat[,2]=time2
  tp=dat
  ooo=order(time2)
  tp=tp[ooo,]
  ## order the data in increasing order of exit time
  ## cat("ELSE loop \n")      

  utimes=unique(tp[,2][tp[,3]==1])  
  ## cat("utimes=",utimes[1:NROW(utimes)], "\n")
  ## this is same as t.evaluate from riskset w/o entry
  utimes=utimes[order(utimes)]
  nutimes=length(utimes)
  ## this is same as nt
  nt=nutimes
  nobs=NROW(tp)
  ## this is same as nc
  nc=nobs
  t.evaluate=c(min(time),utimes)
  newdata=NULL
  tpnew=NULL
  newStatus=NULL
  for(i in 1:nt)
    {
      start=matrix(t.evaluate[i], nc,1)
      finish=matrix(t.evaluate[i+1], nc,1)
      newStatus=ifelse((tp[,2]==t.evaluate[i+1])&(tp[,3]!=0), 1, 0)
      keeps=ifelse((tp[,1]<=t.evaluate[i+1])&(tp[,2]>=t.evaluate[i+1]), 1, 0)
      ncl=ncol(tp)
      tpnew=cbind(start, finish, newStatus, tp[,4:ncl])
      tpnew=tpnew[keeps==1,]
      newdata=rbind(newdata, tpnew)
    }
  rs.tp=newdata
  colnames(rs.tp)=c("start", "finish", "newStatus", colnames(tp)[4:ncl]) 
  return(rs.tp)
}


## 5. weightedKM(Stime, status, wt=NULL, entry=NULL): This function
## estimate S(t) where sampling weights are permitted. 
## *** NOTE: function name changed from my.km to weightedKM as of
## Sep29,2006

## #############################
## THIS FUNCTION IS WORKING A OK
## #############################

weightedKM<-function(Stime, status, wt=NULL, entry=NULL ){
  ##
  ## PURPOSE:  obtain survival function estimate S(t) where
  ##           sampling weights are permitted.
  ##
  ## DATE: Dec 4, 2006, arguments Stime and entry replaced by time and
  ## time2, like Surv() function

  Sstatus <- status
  if( is.null(entry) )  entry <- rep( 0, length(time) )
  
  if( is.null(wt) ) wt <- rep( 1, length(Stime) )
  ##
  dt <- unique( Stime[Sstatus==1] )
  dt <- dt[order(dt)]
  survival <- rep( NA, length(dt) )
  current <- 1.00
  for( j in 1:length(dt) ){
    at.risk <- as.logical( (entry<=dt[j]) & (Stime>=dt[j]) )
    dead <- as.logical( (entry<=dt[j]) & (Stime==dt[j]) & (Sstatus==1) )
    n <- sum(  at.risk*wt  )
    d <- sum(  dead*wt  )
    current <- current*( 1.00 - d/n )
    survival[j] <- current
  }
  out <- list( time=dt, survival=survival )
  out
}

## *** Note: Change Sstatus to status
## *** Note: This is same as survfit() function with appropriate
## options, but it takes longer than survfit(). For example, with PBC1
## data, the times taken are:
## system.time(survfit( Surv(survival.time,survival.status) ~ 1 , conf.type="none"))
## [1] 0.009 0.000 0.009 0.000 0.000
## system.time(weightedKM( Stime=survival.time, status=survival.status ))
## [1] 0.020 0.001 0.021 0.000 0.000


## 6. KM.plot(times, survival, max.T=NULL, lty=NULL, all=T )
## This function plots Kaplan-Meier plot (creates step function).
## changed times to Stime as of Sep29, 2006
## changed the argument list to include a ... argument as discussed in
## Writing R Extension
## function name changed from plot.KM to KM.plot for R issues
KM.plot<-function( Stime, survival, max.T=NULL, lty=NULL, all=TRUE,...){
  ##
  ## PURPOSE:  create Kaplan-Meier plot
  ##
  ##
  times <- Stime
  n <- length( times )
  times<-rep( times, rep(2,n) )
  survival <-rep( survival, rep(2,n) )
  ##
  if( is.null(max.T) ) max.T <- max(times)+1
  if( is.null(lty) ) lty <- 1
  ##
  ##
  ## Offset and Add head and tail
  ##
  ##
  times <- c( 0, times, max.T )
  survival <- c( 1, 1, survival )
  ##
  if( all ){
    plot( times, survival, xlim=c(0,max.T), ylim=c(0,1), lty=lty,
         type="l", xlab="Time", ylab="Survival",... )
  }else{
    lines( times, survival, lty=lty )
  }
  ## end-of-fnx
}
##---------------------------------------------------------------------------

## 7. SchoenSmooth(fit, Stime, status, span=0.40,
##    order=0, shrink=F ) 
##    This function smooth the Schoenfeld residuals using
##    Epanechnikov's optimal kernel. 
## ***Oct 6,2006: argument shrink=F/T is removed

## ################################
## ASK IF THIS NEEDED TO BE CHANGED
## ################################

SchoenSmooth <- function(fit, Stime, status, span=0.40, order=0, entry=NULL)
{
  ##
  shrink=FALSE
  survival.time <- Stime
  survival.status <- status
  phtest <- cox.zph( fit, transform="identity" )
  ##
  ## edited:  18 July 2003
  ## utimes <- unique( survival.time[ survival.status==1 ] )
  utimes <- survival.time[ survival.status==1 ]
  utimes <- utimes[ order(utimes) ]
  
  p <- length( fit$coef )

  if( p==1 )
    {
      bbb <- rep( NA, length(utimes) )
      shrinkage <- rep( NA, length(utimes) )
    }
  else
    {
      bbb <- matrix( NA, length(utimes), p  )
      shrinkage <- matrix( NA, length(utimes), p  )
    }
  ##
  ##
  ## find the bandwidth for each time
  ##
  grid.t <- utimes
  nt <- length(grid.t)
  lambda <- rep(0,nt)
  ##
  ## take span% of the data
  ##
  nd <- round( span*nt )
  ##
  for (i in 1:nt)
    {
      lambda[i] <- sort( abs( grid.t-grid.t[i] ) )[nd]
    }
  ##
  for( j in 1:length(utimes) )
    {
    ##
    target <- grid.t[j]
    ## use Epanechnikov's optimal kernel
    u <- (grid.t - target)/lambda[j]
    iu <- ifelse( abs(u)>1, 0, 1 )
    wt <- iu*0.75*(1-u^2)/lambda[j]
    ##
    www <- wt/sum(wt)
    ##
    if( p==1 )
      {
      if( order==0 )
        {
        beta <- sum( phtest$y * www )
        nnn <- sum( www>0 )
        ddd <- (nnn-1)/nnn
        if( shrink )
          {
          beta <- ddd * beta
        }
      }
      if( order==1 )
        {
        cTime <- phtest$x - target
        fit <- lm( phtest$y ~ cTime, weights=www )
        nnn <- sum( www>0 )
        ddd <- (nnn-1)/nnn
        beta <- fit$coef[1]        
        if( shrink )
          {
          beta <- ddd * fit$coef[1]          
        }
      }
      bbb[j] <- beta
      shrinkage[j] <- ddd
    }
    else
      {
      for( k in 1:p )
        {
        if( order==0 )
          {
          beta <- sum( phtest$y[,k] * www )
          nnn <- sum( www>0 )
          ddd <- (nnn-1)/nnn
          if( shrink )
            {
            beta <- ddd * beta            
          }          
        }
        if( order==1 )
          {
          cTime <- phtest$x - target
          fit <- lm( phtest$y[,k] ~ cTime, weights=www )
          beta <- fit$coef[1]
          nnn <- sum( www>0 )
          ddd <- (nnn-1)/nnn          
          if( shrink ){
            beta <- ddd * beta
          }          
        }
        bbb[j,k] <- beta
        shrinkage[j,k] <- ddd
      }
    }
  }
  ## out <- list( time=utimes, beta=bbb, shrinkage=shrinkage )
  ## DATE: OCT 18, 2006 -- shrinkage REMOVED FROM RETURN LIST
  out <- list( time=utimes, beta=bbb)
  out
}



######################################################



## 8. risksetROC(time, time2=NULL, status, eta, predict.time, method="Cox"
## or "LocalCox" or "Schoenfeld", span=NULL, p=1, order=1,
## window="asymmetric", plot=TRUE, type="l",xlab="FP",ylab="TP",...) 


risksetROC <- function(Stime, entry=NULL, status, marker, predict.time, method="Cox",
                       span=NULL, order=1, window="asymmetric", prop=0.5,
                       plot=TRUE, type="l", xlab="FP", ylab="TP", ...)

  ## DATE: OCT 18, 2006 -- p=1 REMOVED FROM ARGUMENT  LIST
  ## argument eta changed to marker
  ## prop denotes what proportion of length of time-interval to use for llCoxReg()
  {
    ## *** NOTE: eta is M (marker) from the paper
    mChoice <- match(method, c("Cox","LocalCox", "Schoenfeld"))
    if( is.na(mChoice) )
      {
        cat("error in method choice")
        stop(0)
      }
    ## In any case, FP can be obtained from eqn(2) of paper, so we
    ## need Coxweight() to run anyway.
    if(is.null(span)&((mChoice==2)||(mChoice==3)))
      {
        cat("Need span for methods = \"LocalCox\" or \"Schoenfeld\" \n")
        stop(0)
      }
    p=1
    eta=marker
    if(length(entry)==0)
      {
        entry=rep(0,NROW(Stime))
      }
    time=entry
    time2=Stime
    
    if(mChoice==1) ## "Cox"
      {
        fit=coxph(Surv(time, time2, status)~eta)
        new.eta=eta*fit$coefficients
      }
    if(mChoice==2) ## "LocalCox"
      {
        nt=length(unique(time2[status==1]))
        grid.t=time2
        nd=prop
        delta=abs(grid.t-predict.time)/(max(time2[status==1])-min(time2[status==1]))
        keep=delta<=nd

        bfnx.ll=llCoxReg(entry=time[keep==1],Stime=time2[keep==1], status=status[keep==1],
          marker=eta[keep==1], span=span, p=p, window=window)

        gamma.t=bfnx.ll$beta[NROW(bfnx.ll$time[bfnx.ll$time<=predict.time]),]
        new.eta=eta*gamma.t[1]
      }
    if(mChoice==3) ## "Schoenfeld"
      {
        fit=coxph( Surv(time, time2, status) ~ eta )            
        
        bfnx.SS=SchoenSmooth(fit=fit, Stime=time2, status=status, span=span,
          order=order)
        gamma.t=bfnx.SS$beta[NROW(bfnx.SS$time[bfnx.SS$time<=predict.time])]
        new.eta=eta*gamma.t
      }
    out=CoxWeights(marker=new.eta, entry=time, status=status,
      predict.time=predict.time, Stime=time2)
    if(plot==TRUE)
      {
        plot(out$FP, out$TP, type=type, xlab=xlab, ylab=ylab,...)
        abline(c(0,0), c(1,1))
      }
    return(out)
  }


## 9. risksetAUC(time, time2, status, eta, method="Cox"
## or "LocalCox" or "Schoenfeld", span=NULL, p=1, order=1,
## window="asymmetric", plot=TRUE, ...) 

risksetAUC <- function(Stime, entry=NULL, status, marker, method="Cox",
                       span=NULL, order=1, window="asymmetric",
                       tmax, weight="rescale", plot=TRUE, type="l",
                       xlab="Time", ylab="AUC", ...)
##   function(time, time2=NULL, status, marker, method="Cox",
##                        span=NULL, order=1, window="asymmetric",
##                        tmax, weight="rescale", plot=TRUE, type="l",
##                        xlab="Time", ylab="AUC", ...)
  ## DATE: OCT 18, 2006 -- p=1 REMOVED FROM ARGUMENT  LIST
  ## argument eta changed to marker
{
 
  mChoice <- match(method, c("Cox","LocalCox", "Schoenfeld"))
  if( is.na(mChoice) )
    {
      cat("error in method choice")
      stop(0)
    }
  ## In any case, FP can be obtained from eqn(2) of paper, so we
  ## need Coxweight() to run anyway.
  ## identify unique failure times

    if(is.null(span)&((mChoice==2)||(mChoice==3)))
    {
      cat("Need span for methods = \"LocalCox\" or \"Schoenfeld\" \n")
      stop(0)
    }
  p=1
  eta=marker
  if(length(entry)==0)
    {
      time=rep(0, length(Stime))
    }
  else
    {
      time=entry
    }
  Stime=Stime
  time2=Stime
  utimes=unique(Stime[status==1])
  utimes=utimes[order(utimes)]
  
  ## cat(utimes)
  new.eta=NULL

  km.out= weightedKM(Stime=Stime, status=status, entry=time)

  if(mChoice==1) ## "Cox"
    {
      fit=coxph(Surv(time=time, time2=time2, event=status) ~ eta)
      ## #############################################      
      gamma.t=rep(fit$coefficients,NROW(utimes))
      ## check this line --  checked OK :D
      ## NOTE that in risksetROC() we multiply coeff with eta here
      ## but in this function, we multiply coeff with eta later
      ## #############################################
    }

  if(mChoice==2) ## "LocalCox"
    {
      bfnx.ll=llCoxReg(Stime=Stime, status=status, marker=eta, span=span,
        p=p, window=window, entry=entry)
      gamma.t=bfnx.ll$beta[,1]
    }

  if(mChoice==3) ## "Schoenfeld"
    {
      fit=coxph(Surv(time=time, time2=time2, event=status) ~ eta)
      bfnx.SS=SchoenSmooth(fit=fit, Stime=Stime, status=status, span=span,
        order=order)
      gamma.t=bfnx.SS$beta
      gamma.t=gamma.t[!duplicated(gamma.t)]
      ## str(bfnx.SS); str(gamma.t); str(utimes)
    }
  AUC=NULL

  for(i in 1:NROW(gamma.t))
    {
      ## cat("i=", i, "utimes=", utimes[i],"\n")
      new.eta=eta*gamma.t[i]
      out=CoxWeights(marker=new.eta, Stime=Stime, status=status,
    predict.time=utimes[i], entry=entry)
      ## str(out)
      AUC=c(AUC,out$AUC)
    }
  Cindex=IntegrateAUC(AUC=AUC, utimes=utimes, St=km.out$survival,
    tmax=tmax, weight=weight) 

  if(plot==TRUE)
    {
      plot(utimes, AUC, type=type, xlim=c(min(utimes), tmax+1),
    ylim=c(0.4,1.0), xlab=xlab, ylab=ylab,...)
      abline(h=0.5)
    }
  return(out=list(utimes=utimes, St=km.out$survival,
           AUC=AUC, Cindex=Cindex))
  
}
