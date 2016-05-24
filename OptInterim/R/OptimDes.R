OptimDes<-
function(B.init,m.init,alpha,beta,param,x,target=c("EDA","ETSL","ES"),
         sf=c("futility","OF","Pocock"),num.arm,r=0.5,num.stage=2,
         pause=0,
         control=OptimDesControl(),...)
{

  ##set the seed for reproducibility of function optim
  if(exists('.Random.seed')){
    oldseed<-.Random.seed
    set.seed(12357)
  }else{
    set.seed(12357)
    oldseed<-.Random.seed
  }
  on.exit(.Random.seed<-oldseed)


  target <- match.arg(target)
  sf <- match.arg(sf)

  if(!num.arm%in%c(1,2))
  {
    stop("the number of treatment arms must be either one or two")
  }   

  shape0 <- param[1]
  scale0 <- param[2]
  shape1 <- param[3]
  scale1 <- param[4]

  trace <- control$trace
  tol <- control$tol
  startC1L <- control$startC1L
  n.int<-control$n.int[num.arm]
  aboveMin<- control$aboveMin[num.arm]

  if(identical(num.arm,1) && !identical(sf,"futility"))
  {
    stop("Only futility interim analysis for single-arm trial")
  }
  
  if(!num.stage%in%c(2,3))
  {
    stop("The optimal design allows two or three stages")
  }

  if(identical(num.arm,1) && !identical(r,0.5))
  {
    warning("Randomization ratio not equal to 0.5 for a single-arm study")
  }
  
  if(identical(num.arm,2) && isTRUE(r*(1-r)<=0))
  {
    stop("invalid randomization ratio r")
  }
  
    
  if(!identical(length(B.init),length(m.init)))
  {
    stop("Projected patient times and numbers should be of equal length")
  }

  if(isTRUE(x>=max(B.init)))
  {
    stop("The survival time of interest is beyond the projected accrual time")
  }

  if(isTRUE(pause>x))stop("Pauses in accrual should not exceed the target time")


  if(!identical(target,"EDA") && !identical(target,"ETSL") && 
      !identical(target,"ES"))
    stop("unknown target value")
  
  if(!identical(trace,FALSE) && !identical(trace,TRUE))
    stop("trace should be a logical value.")
    
  ## Step 1, input n, rho1
  b <- length(B.init)

  p0<-s(shape0,scale0,x)
  p1<-s(shape1,scale1,x)
  if(isTRUE(p1<=p0))
  {
    stop("the alternative survival rate should be greater than the null survival rate")
  }

  # single stage sample size, duration of accrual and study length
  fix.d <- FixDes(B.init,m.init,alpha,beta,param,x,num.arm,r)
  n0 <- fix.d$n0
  da <- fix.d$DA
  sl <- fix.d$SL

  ## compute fixed sample size based on exact distributions for
  ## use in optional proportional adjustment of sample sizes and times
  n0E<- fix.d$n0E
  DAE<-fix.d$DAE
  SLE<-fix.d$SLE

  if(isTRUE(n0>floor(sum(m.init))))
  {
    stop(paste("projected patient sample size is below the minimum requirement of ", 
                n0))
  }
  if(isTRUE(x>=da && target!='ETSL'))stop(paste('Single stage design is optimal because all patients can be',
  'accrued before the event time'))

  # max sample sizes searched
  nsearch<-seq(n0,sum(m.init),n.int)
  nvec<-length(nsearch)

  # store  information
  nrho.res <- matrix(rep(NA,nvec*3),nrow=nvec,ncol=3)
  nrho.res[,1] <- nsearch
  n.res <- nsearch
  nvec<-length(n.res)
  # set up vector to hold the minimum EDA, ES or ETSL for each n
  min.res <- rep(NA,nvec)
  t1.res <- rep(NA,nvec)
  t2.res <- rep(NA,nvec)
  C1.res <- rep(NA,nvec)
  C1U.res <- rep(NA,nvec)
  C2.res <- rep(NA,nvec)
  C2U.res <- rep(NA,nvec)
  C3.res <- rep(NA,nvec)
  MTSL.res <- rep(NA,nvec)
  MDA.res <- rep(NA,nvec)
  ETSL.res <- rep(NA,nvec)
  EDA.res <- rep(NA,nvec)
  ES.res <- rep(NA,nvec)
  pstopNull.res <- rep(NA,nvec)
  pstopAlt.res <- rep(NA,nvec)
  n1.res <- rep(NA,nvec)
  n2.res <- rep(NA,nvec)

  if(identical(num.stage,2))
  {
  u.res <-  matrix(rep(NA,2*nvec),ncol=2)       
  se.res <- matrix(rep(NA,4*nvec),ncol=4)
  }
  if(identical(num.stage,3))
  {
  u.res <-  matrix(rep(NA,3*nvec),ncol=3)       
  se.res <- matrix(rep(NA,6*nvec),ncol=6)
  }  
  
  currentMin<-Inf
  i<-0
  for(n in nsearch)
  {
    i<-i+1
    if(num.stage==2){
        optout <- optimize(f=f.DesWrap,interval=c(0,1),tol=tol,B.init=B.init,
                m.init=m.init,alpha=alpha,beta=beta,param=param,x=x,n=n,
                num.arm=num.arm,r=r,num.stage=num.stage,target=target,sf=sf,
                pause=pause)          
       
        min.res[i] <- optout$objective
        nrho.res[i,2] <- optout$minimum
        if(trace){
          cat("n=",nrho.res[i,1],"optimal rho=",nrho.res[i,2],
              target,"=",min.res[i],"\n")
          flush.console()
        }

        if(optout$objective==1e10) {
           warning(paste('Optimizaton failed for n= ',nrho.res[i,1]))
           next
        }

    }else{ ### 3 stage design
       ### parameters maximized are (transformed) rho1, rho2, and c1l
       rho1<- 0.5                
       rho2<- 0.75              
       logitrho1<-qlogis(rho1)
       logitdif<-qlogis(rho2)-logitrho1
       loglogitdif<-log(logitdif)
       startpar<-c(logitrho1,loglogitdif,startC1L)
       optout <- optim(par=startpar,fn=f.DesWrap,method="BFGS",
                 B.init=B.init,m.init=m.init,alpha=alpha,beta=beta,
                 param=param,x=x,
                 n=n,num.arm=num.arm,r=r,num.stage=num.stage,
                 target=target,sf=sf,pause=pause)

       min.res[i] <- optout$value
       rho1<-plogis(optout$par[1])
       rho2<-plogis( optout$par[1]+exp(optout$par[2]) )
       nrho.res[i,2] <- rho1
       nrho.res[i,3] <- rho2
       C1.res[i] <- optout$par[3]
       if(trace){
          cat("n=",nrho.res[i,1],"optimal rho=",nrho.res[i,-1],
              target,"=",min.res[i]," convergence=",optout$convergence,"\n")
          flush.console()
       }
       
       if(optout$convergence!=0 || optout$value==1e10){
         min.res[i] <- Inf
         nrho.res[i,2] <- NA
         nrho.res[i,3] <- NA
         C1.res[i] <- NA
         warning(paste('Optimizaton failed for n= ',nrho.res[i,1]))
         next
       }     
    }

    ### re-run final converged values to obtain full output
    if(num.stage==2) {optpar <- nrho.res[i,2]
    }else if(num.stage==3) optpar <-c(nrho.res[i,2:3],C1.res[i])

    f.out <- f.Des(optpar,B.init,m.init,alpha,beta,param,x,n,
                   num.arm,r,num.stage,target,sf,pause=pause)
    if(f.out$min==1e10){ 
        min.res[i] <- Inf
        next
    } 

    MDA.res[i] <- f.out$mda
    MTSL.res[i] <- f.out$mda+x
    se.res[i,] <- f.out$se
    u.res[i,] <- f.out$u
    EDA.res[i]<-f.out$EDA
    ETSL.res[i] <- f.out$ETSL
    ES.res[i] <- f.out$ES
    pstopNull.res[i]<-f.out$pstopNull
    pstopAlt.res[i]<-f.out$pstopAlt
    n1.res[i] <- ceiling(f.out$n1)
    t1.res[i] <- f.out$t1
    C1U.res[i] <- f.out$C1U

    if(num.stage==2){
      C1.res[i] <- f.out$C1
      C2.res[i] <- f.out$C2
    }else if(num.stage==3){ 
      t2.res[i] <- f.out$t2
      n2.res[i] <- ceiling(f.out$n2)
      C2.res[i] <- f.out$C2L
      C2U.res[i] <- f.out$C2U
      C3.res[i] <- f.out$C3
    }

    if(min.res[i]<currentMin){currentMin<-min.res[i]
    }else if(min.res[i]>aboveMin*currentMin)break
  }

  ### return an error message if optimization fails for all n evaluated
  if(!is.finite(currentMin))stop('An optimal design could not be computed')

  ## truncate result vectors at last sample size evaluated
  nrho.res <- nrho.res[1:i,]
  min.res <-  min.res[1:i]
  n.res <-    n.res[1:i]
  t1.res <-   t1.res[1:i]
  t2.res <- t2.res[1:i]
  C1.res <-   C1.res[1:i]
  C1U.res <- C1U.res[1:i]
  C2.res <-   C2.res[1:i]
  C2U.res <- C2U.res[1:i]
  C3.res <- C3.res[1:i]
  MTSL.res <- MTSL.res[1:i]
  MDA.res <-  MDA.res[1:i]
  ETSL.res <- ETSL.res[1:i]
  EDA.res <-  EDA.res[1:i]
  ES.res <-   ES.res[1:i]
  pstopNull.res <-   pstopNull.res[1:i]
  pstopAlt.res <-   pstopAlt.res[1:i]
  n1.res <-   n1.res[1:i]
  n2.res <- n2.res[1:i]
  se.res <-   se.res[1:i,]
  u.res <-  u.res[1:i,]
  

  ## select the optimal n
  order.min <- order(min.res)[1]
  n.last <- nrho.res[order.min,1]
  outcome <- min.res[order.min]
  EDA<-EDA.res[order.min]
  ETSL<-ETSL.res[order.min]
  ES <- ES.res[order.min]
  pstopNull <- pstopNull.res[order.min]
  pstopAlt  <- pstopAlt.res[order.min]
  mda <- MDA.res[order.min]
  t1.last <- t1.res[order.min]
  C1.last <- C1.res[order.min]
  C1U.last <- C1U.res[order.min]
  C2.last <- C2.res[order.min]
  n1 <- n1.res[order.min]
  se <- se.res[order.min,]
  u  <- u.res[order.min,]

  if(t1.last>=mda){
     pamt<-min(pause,mda+x-t1.last)
  }else pamt<-pause
  tinterim<-t1.last+pamt
  # truncC has shifted accrual distribution
  EW<-truncC(B.init,m.init,n.last,x,mda,tinterim,
       pause,pt1=t1.last) #potential exposure at t1
  EW<-n1*EW
  if(identical(num.arm,1) && identical(num.stage,2)){
      if(isTRUE(ceiling((n0E/n0)*n.last)>floor(sum(m.init)))){ 
        warning("The maximum patient accrual is below the minimum requirement for the normal approximation adjustment")
        EWadj<-NA
      }else{
          CMadjfac<-n0E/n0
          mdaAdj<-ifelse(t1.last<mda,CMadjfac*(mda-pause)+pause,CMadjfac*mda)
          tadj<-CMadjfac*t1.last
          tinterim<-tadj+pamt
          EWadj<-truncC(B.init,m.init,ceiling(CMadjfac*n.last),
          x,mdaAdj,tinterim, pause,pt1=tadj)
          EWadj<-n1*(n0E/n0)*EWadj
          EW<-c(EW,EWadj)
      }
  }

  ## compare multi-stage to single-stage
  singleEff<-n0
  if(identical(target,'EDA'))singleEff<-da
  if(identical(target,'ETSL'))singleEff<-sl
  if(isTRUE(singleEff<=currentMin))warning('The single stage design was more efficient than the multi-stage design')

  if(identical(num.stage,2)){

      res <-structure(list(target=target,sf=sf,test=c(alpha=alpha,beta=beta,
            param=param,x=x), 
            design=c(num.arm=num.arm,r=r,num.stage=num.stage,pause=pause),
            accrual=list(B.init=B.init,m.init=m.init),
            result=c(EDA=EDA,ETSL=ETSL,ES=ES,pstopNull=pstopNull,pstopAlt=pstopAlt),
            n=c(n1=n1,n.last=n.last),
            stageTime=c(t1=t1.last,MTSL=mda+x,MDA=mda),boundary=c(C1L=C1.last,
            C1U=C1U.last,C2=C2.last),se=se,u=u,exposure=EW,
            all.info=data.frame(n=n.res,t1=t1.res,C1L=C1.res,C1U=C1U.res,C2=C2.res,
            MDA=MDA.res,MTSL=MTSL.res,EDA=EDA.res,ETSL=ETSL.res,ES=ES.res),
            single.stageTime=c(n0=n0,DA=da,SL=da+x,n0E=n0E,DAE=DAE,SLE=SLE)),
            class="OptimDes")
  }else if(identical(num.stage,3)){
      t2.last <- t2.res[order.min]
      C2U.last <- C2U.res[order.min]
      C3.last <- C3.res[order.min]
      n2 <- n2.res[order.min]

      if(t2.last>=mda){
         pamt2<-min(pause,mda+x-t2.last)
      }else pamt2<-pause
      tinterim<-t2.last+pamt2
      EW2<-truncC(B.init,m.init,n.last,x,mda,tinterim,pause=pause,
                  pt1=t1.last,pt2=t2.last) #potential exposure at t2
      EW<-c(EW,n2*EW2)   # add exposure at second interim analysis

      res <-structure(list(target=target,sf=sf,test=c(alpha=alpha,beta=beta,
            param=param,x=x), 
            design=c(num.arm=num.arm,r=r,num.stage=num.stage,pause=pause),
            accrual=list(B.init=B.init,m.init=m.init),
            result=c(EDA=EDA,ETSL=ETSL,ES=ES,pstopNull=pstopNull,pstopAlt=pstopAlt),
            n=c(n1=n1,n2=n2,n.last=n.last),
            stageTime=c(t1=t1.last,t2=t2.last,MTSL=mda+x,MDA=mda),
            boundary=c(C1L=C1.last,
            C1U=C1U.last,C2L=C2.last,C2U=C2U.last,C3=C3.last),se=se,
            u=u,exposure=EW,
            all.info=data.frame(n=n.res,t1=t1.res,t2=t2.res,C1L=C1.res,
            C1U=C1U.res,C2L=C2.res,
            C2U=C2U.res,C3=C3.res,MDA=MDA.res,MTSL=MTSL.res,EDA=EDA.res,
            ETSL=ETSL.res,ES=ES.res),
            single.stageTime=c(n0=n0,DA=da,SL=da+x,n0E=n0E,DAE=DAE,SLE=SLE)),
            class="OptimDes")
  }  
  return(res)   
}

##############################################################################################  

  ### functions for the weibull distribution
  ### h and s use R definitions of weibull parameters

  # hazard function of Weibull distribution 
  h <- function(shape,scale,x)
  {
    (shape/scale)*(x/scale)^(shape-1)
  }

  # survival function of Weibull distribution
  s <- function(shape,scale,x)
  {
    exp(-(x/scale)^shape)
  }

  # cummulative hazard function of weibull distribution
  lambda <- function(shape,scale,x)
  {
    (x/scale)^shape
  }
  
##############################################################################################

  ### accrual distribution function

  # function for Fn(t)
  fnt <- function(B.init,m.init,n,t,mda,l)
  {
    B <- B.init[1:l]   # "L" not one
    B[l] <- mda
    M <- m.init[1:l]
    M[l] <- ifelse(l>1,n-sum(m.init[1:(l-1)]),n)
    ll <- length(t)
    m1 <- matrix(rep(B,ll),nrow=l,byrow=FALSE)
    m2 <- matrix(0,nrow=l+1,ncol=ll)
    m2[1:l,] <- m1
    m2[l+1,] <- t
    # the interval that t falls in    
    count <- apply(m2,2,indX)
    M.sum <- cumsum(M)
    count.1 <- ifelse(count==1,1,count-1) # avoid error message in the following ifelse function
    pr <-ifelse(count==1,(M[1]/n)*(t/B[1]),ifelse(count<=l,M.sum[count.1]/n+
         (M[count]/n)*((t-B[count.1])/(B[count]-B[count.1])),1))
    pr<-ifelse(t>=mda,1,pr)
    return(pr)
      
  }    


  ### function to find smallest index for entry exceeding an input

  indX<-function(x) {
    rank(x,ties.method="min")[length(x)]
  }

  ### 
    compMDA<-function(B.init,m.init,n){
       l<-indX(c(cumsum(m.init),n))
       mda <- ifelse(l>1,B.init[l-1]+(B.init[l]-B.init[l-1])*(n-
            sum(m.init[1:(l-1)]))/m.init[l],B.init[l]*(n/m.init[l]))
       return(c(mda,l))
    }
    compTime<-function(B.init,m.init,ti){
       l1<-indX(c(B.init,ti))
       n1 <- ifelse(l1>1,sum(m.init[1:(l1-1)])+ m.init[l1]*(ti-
              B.init[l1-1])/(B.init[l1]-B.init[l1-1]),
              m.init[1]*ti/B.init[1])    
       return(n1)
    }

  ### compute the single stage design sample size using exact binomial for one-arm study
  ## one arm
  single.exact<-function(n0,alpha,beta,p0,p1){
     if(p0>=p1)stop("The null event-free rate exceeds the alternative rate")
     fixN<-F
     n0E<-floor(n0/2)
     while(!fixN){

        bx<-qbinom(alpha,n0E,p0,lower.tail=F)
        if( pbinom(bx,n0E,p1)<=beta){
           fixN<-T
        }else n0E<-n0E+1

     }
  return(n0E)
  }

truncC<-function(B.init,m.init,n,x,mdaP,ti,pause=0,pt1=NULL,pt2=NULL){
  ### Form distribution of Y given Y<=ti
  ### then return the expectation of  the conditional
  ### exposure time

  
  mmod<-m.init
  Bmod<-B.init

  if((!is.null(pt1)) && (pt1<mdaP)){
    outshift<-shiftDist(B.init,m.init,x,pause,pt1)
    mmod<-outshift$mmod
    Bmod<-outshift$Bmod
  } 
  if((!is.null(pt2)) && (pt2<mdaP)){
    outshift<-shiftDist(Bmod,mmod,x,pause,pt2)
    mmod<-outshift$mmod
    Bmod<-outshift$Bmod
  }
  
  ### compute conditional density
  k<-1
  while(Bmod[k]<ti)k<-k+1

  m<-mmod[1:k]
  B<-c(0,Bmod[1:k])  
  m[k]<-mmod[k]*(ti-B[k])/(B[k+1]-B[k])
  B[k+1]<-ti

  Pi<-m/sum(m)
  dB<-diff(B)
  Di<-ifelse(dB>0,Pi/diff(B),0)

  ### compute conditional expectation
  EW<-truncM(B,Di,x,ti)
  return(EW)
}

shiftDist<-function(B,m,x,pause,pt){
  ## shift distribution for pause in accrual

  ## no shift if pt after end of accrual
  if(pt>=max(B))return(list(mmod=m,Bmod=B))

  k<-1
  while(B[k]<pt)k<-k+1
 
  bm<-B[1:k]
  mm<-m[1:k]

  Bmod<-c(0,B)

  ### truncate inteval containing pt
  bm[k]<-pt
  mm[k]<-m[k]*(pt-Bmod[k])/(Bmod[k+1]-Bmod[k])
  ### add pause without accrual
  bm<-c(bm,pt+pause)
  mm<-c(mm,0)
  ### if not on boundary, add remainder of the interval
  if(pt<B[k]){
      mm<-c(mm,m[k]*(Bmod[k+1]-pt)/(Bmod[k+1]-Bmod[k]))
      bm<-c(bm,bm[k+1]+(Bmod[k+1]-pt))
  }
  ### add in remainder of the accrual distribution
  nb<-length(B)
  if(k<nb){
    mm<-c(mm,m[(k+1):nb])
    bm<-c(bm,B[(k+1):nb]+pause)
  }
  return(list(mmod=mm,Bmod=bm))
}

truncM<-function(B,Di,x,t1){
  ### Compute the expectation of exposure over the conditional distribution
  trM<-0
  b<-length(B)
  for (i in 2:b){

     if(t1-x>=B[i]){
         trM<-trM + x*Di[i-1]*(B[i]-B[i-1])
     }else if(t1-x>=B[i-1]){
         tmp<-x*(t1-x-B[i-1]) + t1*(B[i]-t1+x) - 0.5*(B[i]^2-(t1-x)^2)
         trM<-trM+Di[i-1]*tmp
     }else{
         tmp<- t1*(B[i]-B[i-1]) - 0.5*(B[i]^2-B[i-1]^2)
         trM<- trM+Di[i-1]*tmp
     }

  }
  return(trM)
}


  # function for equation 1
  sig2 <- function(B.init,m.init,shape,scale,x,n,t,mda,l,
                   pause=0,pt1=NULL,pt2=NULL)
  ### up to two pauses can be added, accrual and exposure
  ### adjusted accordingly.  pt1 and pt2 are the times
  ### when the accrual pauses begin
  {
      sig2.int <- integrate(f.int,lower=0,upper=x,
                  B.init=B.init,m.init=m.init,shape=shape,
                  scale=scale,x=x,n=n,t=t,mda=mda,l=l,
                  pause=pause,pt1=pt1,pt2=pt2
                  )
      if(sig2.int$message!="OK")
        stop("integration for sig2 cannot be completed")
      else
        sig2.int$value
  }
### function to be integrated in function sig2
f.int <- function(a,B.init,m.init,shape,scale,x,n,t,mda,l,
                   pause,pt1,pt2){
  ta<-t-a
  if(isTRUE(all.equal(pause,0)) || isTRUE(pt1>mda)){ ###no prior pauses
    pac<-fnt(B.init,m.init,n,ta,mda,l)
  }else{ 
    pac<-numeric(length(ta))
    tc1<-(ta<=pt1)
    pac[tc1]<-fnt(B.init,m.init,n,ta[tc1],mda,l)

    tc2<-(!tc1) & (ta<pt1+pause) 
    pac[tc2]<-fnt(B.init,m.init,n,pt1,mda,l)

    tc3<-(!tc1) & (!tc2)
    pac[tc3]<-fnt(B.init,m.init,n,ta[tc3]-pause,mda,l)

    if(!is.null(pt2)){
        if((pt2>pt1) && (pt2<=mda+pause)){  ## 2nd interim before accrual end
            tc4<-(ta>pt2) & (ta<=pt2+pause)
            pac[tc4]<-fnt(B.init,m.init,n,pt2-pause,mda,l)
            tc5<-(ta>pt2+pause)
            pac[tc5]<-fnt(B.init,m.init,n,ta[tc5]-2*pause,mda,l)
        }
    }
  }
  h(shape,scale,a)/(s(shape,scale,a)*pac)
}

###########################################################################
#### functions for computing t1 and t2
f51G <- function(t,B.init,m.init,shape1,scale1,x,n,mda,l,
                   pause,sig11){
       ### evaluate info after pause
       if(t>=mda){
          pamt<-min(pause,mda+x-t)
       }else pamt<-pause

       sig2(B.init,m.init,shape1,scale1,x,n,t+pamt,
       mda,l,pause=pause,pt1=t)-sig11^2
}

f52G <-  function(t,B.init,m.init,shape1,scale1,shape0,scale0,
                   x,n,mda,l,r,pause,sig1,sig0,lam1,lam0,rho1) {
       if(t>=mda){
          pamt<-min(pause,mda+x-t)
       }else pamt<-pause

      sqrt(sig2(B.init,m.init,shape0,scale0,x,n,t+pamt,
      mda,l,pause=pause,pt1=t)/((1-r)*(lam0^2))+
      sig2(B.init,m.init,shape1,scale1,x,n,t+pamt,
      mda,l,pause=pause,pt1=t)/(r*(lam1^2)))-
      (1/rho1)*sqrt((sig0^2)/((1-r)*(lam0^2))+(sig1^2)/(r*(lam1^2)))
}

f51G2 <- function(t,B.init,m.init,shape1,scale1,x,n,mda,l,
                   pause,t1,sig21){
       ### evaluate info after pause
       if((t1>=mda) || (t>=mda+pause)){
          pamt<-min(pause,mda+x-t)
       }else pamt<-pause

       sig2(B.init,m.init,shape1,scale1,x,n,t+pamt,
       mda,l,pause=pause,pt1=t1,pt2=t)-sig21^2
}

f52G2 <-  function(t,B.init,m.init,shape1,scale1,shape0,scale0,
                   x,n,mda,l,r,pause,t1,sig31,sig30,lam1,lam0,rho21) {

       if((t1>=mda) || (t>mda+pause)){
          pamt<-min(pause,mda+x-t)
       }else pamt<-pause

      sqrt(sig2(B.init,m.init,shape0,scale0,x,n,t+pamt,
      mda,l,pause=pause,pt1=t1,pt2=t)/((1-r)*(lam0^2))+
      sig2(B.init,m.init,shape1,scale1,x,n,t+pamt,
      mda,l,pause=pause,pt1=t1,pt2=t)/(r*(lam1^2)))-
      (1/rho21)*sqrt((sig30^2)/((1-r)*(lam0^2))+(sig31^2)/(r*(lam1^2)))
}

###########################################################################
##### functions for solving for alpha and power
        ### normal integrands for C2L 

        fC2U <- function(C2U,C1U,sigma0,alphaI1,alphaI2){
                 pmvnorm(lower=c(-Inf,C2U),upper=c(C1U,Inf),
                 sigma=sigma0[1:2,1:2])-(alphaI2-alphaI1)
        }

        fC3U <- function(C3U,C1U,C2U,sigma0,alpha,alphaI2){
                 pmvnorm(lower=c(-Inf,-Inf,C3U),upper=c(C1U,C2U,Inf),
                 sigma=sigma0)-(alpha-alphaI2)
        }


        powC1<-function(C1,C1U,C2,sigma1,beta,u1,u2){  # C1 is C1L in the two-arm
                1-pnorm(C1U-u1)+pmvnorm(lower=c((C1-u1),(C2-u2)),
                upper=c((C1U-u1),Inf),sigma=sigma1)-(1-beta)
            }

        powC2<-function(C2,C1U,sigma0,alpha,alphaI){
            1-pnorm(C2)- pmvnorm(lower=c(C1U,C2),
                  upper=c(Inf,Inf),sigma=sigma0)- (alpha-alphaI)
        }

        powC2L<-function(C2L,C1L,C1U,C2U,C3,sigma1,beta,u1,u2,u3){  
                1-pnorm(C1U-u1)+pmvnorm(lower=c((C1L-u1),(C2U-u2)),
                upper=c((C1U-u1),Inf),sigma=sigma1[1:2,1:2])+
                pmvnorm(lower=c((C1L-u1),(C2L-u2),C3-u3),
                upper=c((C1U-u1),(C2U-u2),Inf),sigma=sigma1) -
                (1-beta)
            }



##############################################################################################
  f.Des <- function(optpar,B.init,m.init,alpha,beta,param,x,n,num.arm,r,num.stage,
                    target=c("EDA","ETSL","ES"),sf=c("futility","OF","Pocock"),
                    pause=0)
  {
  # Two-stage: for input n and rho1, derive C1L, C1U and C2 from 
  # the type I&II error constraints
  # Three-stage: for input n, rho11, rho21, C1L, derive C1U and C2L, C2U, C3 
  # from the type I&II error constraints

    target <- match.arg(target)
    sf <- match.arg(sf)

    shape0 <- param[1]
    scale0 <- param[2]
    shape1 <- param[3]
    scale1 <- param[4]    

    lam0<-lambda(shape0,scale0,x)
    lam1<-lambda(shape1,scale1,x)

    smallC<-1e-5       ### small offset 

    if(num.stage==2 && length(optpar)>1)
      stop("For two-stage designs, optpar = rho1")
      
    if(num.stage==3 && length(optpar)!=3)
      stop("For three-stage designs, optpar = c(rho11,rho21,C1L)")

    # optpar is rho1 for two-stage, and c(rho11,rho21,C1L) for three-stage
    if(max(optpar[1:(num.stage-1)]>=1-smallC)) return(list(min=1e10)) ### avoid numerical boundary problems

    ## Step 2, Calculate MDA by n and the accrual function 
    ##         without adjustment for any pauses
    mdaout<-compMDA(B.init,m.init,n)
    mda<-mdaout[1]  ### mda withtout pauses
    l<-mdaout[2]

    if(num.stage==2){
        rho1 <- optpar
        ## Step 3
        ### accrual pauses can be ignored at planned completion
        ### because complete data is available for all patients
        sig21 <- sqrt(sig2(B.init,m.init,shape1,scale1,x,n,mda+x,mda,l))
        sig20 <- sqrt(sig2(B.init,m.init,shape0,scale0,x,n,mda+x,mda,l))
     
        ## Steps 4 and 5
        if(num.arm==1){
            sig11 <- sig21/rho1

            ### find design associated with rho  
            ### function supplied to root finder for step 5 
            
            ### check if uniroot has a solution (largest possible function value should be positive)
            tmp <- f51G(t=x+smallC-pause,B.init=B.init,m.init=m.init,shape1=shape1,
                   scale1=scale1,x=x,n=n,mda=mda,l=l,pause=pause,sig11=sig11)
            if(tmp<0)  return(list(min=1e10))
            t1.root <- try(uniroot(f51G,interval=c(x+smallC-pause,mda+x),
                   B.init=B.init,m.init=m.init,shape1=shape1,
                   scale1=scale1,x=x,n=n,mda=mda,l=l,pause=pause,sig11=sig11
                   ),silent=TRUE)
        }else{
            ## find design associated with rho1
            ### check if uniroot has a solution (largest possible function value should be positive)
            tmp <- f52G(t=x+smallC-pause,B.init=B.init,m.init=m.init,shape1=shape1,
                   scale1=scale1,shape0,scale0,x=x,n=n,mda=mda,
                   l=l,r=r,pause=pause,sig1=sig21,sig0=sig20,
                   lam1=lam1,lam0=lam0,rho1=rho1
                   )
            if(tmp<0)  return(list(min=1e10))            
            t1.root <- try(uniroot(f52G,interval=c(x+smallC-pause,mda+x),
                   B.init=B.init,m.init=m.init,shape1=shape1,
                   scale1=scale1,shape0,scale0,x=x,n=n,mda=mda,
                   l=l,r=r,pause=pause,sig1=sig21,sig0=sig20,
                   lam1=lam1,lam0=lam0,rho1=rho1
                   ),silent=TRUE)
        } 
        if(class(t1.root)=="try-error"){
            stop("the algorithm of uniroot() does not converge in 'maxiter' steps")
        }else{ 
           t1 <- t1.root$root 
        } 

        ### mda adjusted for pause
        if(t1<mda)mdaP<-mda+pause else mdaP<-mda
        ### pause truncated for end of study
        if(t1>=mda){
           pamt<-min(pause,mda+x-t1)
        }else pamt<-pause


        sig10 <- sqrt(sig2(B.init,m.init,shape0,scale0,x,n,t1+pamt,
                 mda,l,pause=pause,pt1=t1))
        if(num.arm==2)sig11 <- sqrt(sig2(B.init,m.init,shape1,scale1,x,
                 n,t1+pamt,mda,l,pause=pause,pt1=t1))


        ## Step 7
        rho0 <- sig20/sig10

        ## Step 8
        if(num.arm==1){
            u2 <- sqrt(n)*(log(lam0)-log(lam1))*lam1/sig21
            u1 <- rho1*u2
            u <- c(u1,u2)
        }else{
            v10 <- (sig10^2)/((1-r)*(lam0^2))
            v11 <- (sig11^2)/(r*(lam1^2))
            v20 <- (sig20^2)/((1-r)*(lam0^2))
            v21 <- (sig21^2)/(r*(lam1^2))
            u1 <- sqrt(n)*(log(lam0)-log(lam1))/sqrt(v10+v11)
            u2 <- sqrt(n)*(log(lam0)-log(lam1))/sqrt(v20+v21) 
            
            u <- c(u1,u2)
        }
        # the number of patients accrued at the time of interim analysis
        tinterim<-ifelse(t1>=mda,mda,t1)  
        n1<-compTime(B.init,m.init,tinterim)

        ## Step 9
        sigma0 <- diag(2)
        sigma1 <- diag(2)
        sigma0[1,2] <- rho0
        sigma0[2,1] <- sigma0[1,2]
        sigma1[1,2] <- rho1
        sigma1[2,1] <- sigma1[1,2]

        alphaI<-0
        C1U <- Inf   
        if(sf != 'futility'){
            alphaI<-spfun(rho0,alpha,sf)
            C1U <- qnorm(1-alphaI)
        }

        C1low<- -3.5      
        C1up<-min(C1U,qnorm(1-alpha))

        if(sf == 'futility'){
            C2<- qnorm(1-alpha)
        }else{
            C2up <- qnorm(1-(alpha-alphaI))
            C2root <- try(uniroot(powC2,c(-3.5,C2up),
                          C1U=C1U,sigma0=sigma0,alpha=alpha,alphaI=alphaI
            ),silent=TRUE)
            if(class(C2root)=="try-error")return(list(min=1e10)) 
            if(suppressWarnings(warning()!="")){
               return(list(min=1e10))  #when there is no solution
            }else{
               C2 <- C2root$root
            }        
        }


        if(sign(powC1(C1low,C1U,C2,sigma1,beta,u1,u2))==
           sign(powC1(C1up,C1U,C2,sigma1,beta,u1,u2))){
           return(list(min=1e10))  #when there is no solution
        }else{
           ###C1 is C1L in the 3-stage trial 
           C1root <- try(uniroot(powC1,c(C1low,C1up),
           C1U=C1U,C2=C2,sigma1=sigma1,beta=beta,u1=u1,u2=u2
           ),silent=TRUE)  
           if(class(C1root)=="try-error")return(list(min=1e10)) 
        }
        if(suppressWarnings(warning()!="")){
           return(list(min=1e10))  #when there is no solution
        }else if(C1root$root>C1U){
           return(list(min=1e10))  #invalid parameter region
        }else{
           C1 <- C1root$root
        }

        ###########################
        pstopNull<-pnorm(C1)+(1-pnorm(C1U))
        pstopAlt<-pnorm(C1-u1)+(1-pnorm(C1U-u1))

        EDA<- ifelse(t1<mdaP,pstopNull*t1+(1-pstopNull)*mdaP,mdaP)
        ETSL<- pstopNull*(t1+pamt)+(1-pstopNull)*(mdaP+x)
        ES<- n1+(1-pstopNull)*(n-n1)
    ## end of two-stage algorithm

    }else if(num.stage==3){
        rho11 <- optpar[1]
        rho21 <- optpar[2]
        C1L <- optpar[3]
        
        if(rho11>=rho21-smallC) return(list(min=1e10))
       
        ## Step 3
        ### accrual pauses can be ignored at planned completion
        ### because complete data is available for all patients
        sig31 <- sqrt(sig2(B.init,m.init,shape1,scale1,x,n,mda+x,mda,l))
        sig30 <- sqrt(sig2(B.init,m.init,shape0,scale0,x,n,mda+x,mda,l))
     
        ## Steps 4 and 5

       if(num.arm==1){
            sig11 <- sig31/rho11

            ### find design associated with rho  
            ### function supplied to root finder for step 5 
            ### check if uniroot has a solution (largest possible function value should be positive)
            tmp <- f51G(t=x+smallC-pause,B.init=B.init,m.init=m.init,shape1=shape1,
                   scale1=scale1,x=x,n=n,mda=mda,l=l,pause=pause,sig11=sig11
                   )
            if(tmp<0)  return(list(min=1e10))
            t1.root <- try(uniroot(f51G,interval=c(x+smallC-pause,mda+x),
                   B.init=B.init,m.init=m.init,shape1=shape1,
                   scale1=scale1,x=x,n=n,mda=mda,l=l,pause=pause,sig11=sig11
                   ),silent=TRUE)
        }else{
            ## find design associated with rho11, rho21
            ### check if uniroot has a solution (largest possible function value should be positive)
            tmp <- f52G(t=x+smallC-pause,B.init=B.init,m.init=m.init,shape1=shape1,
                   scale1=scale1,shape0,scale0,x=x,n=n,mda=mda,
                   l=l,r=r,pause=pause,sig1=sig31,sig0=sig30,
                   lam1=lam1,lam0=lam0,rho1=rho11
                   )
            if(tmp<0)  return(list(min=1e10))        
            t1.root <- try(uniroot(f52G,interval=c(x+smallC-pause,mda+x),
                   B.init=B.init,m.init=m.init,shape1=shape1,
                   scale1=scale1,shape0,scale0,x=x,n=n,mda=mda,
                   l=l,r=r,pause=pause,sig1=sig31,sig0=sig30,
                   lam1=lam1,lam0=lam0,rho1=rho11
                   ),silent=TRUE)
        }
        if(class(t1.root)=="try-error")return(list(min=1e10)) 
        if(suppressWarnings(warning()!="")){
          stop("uniroot to compute t1 failed in 3-stage algorithm")
        }else{ 
             t1 <- t1.root$root 
        } 

        ### mda adjusted for pause at first interim analysis
        if(t1<mda)mdaP<-mda+pause else mdaP<-mda
        ### pause truncated for end of study
        if(t1>=mda){
           pamt<-min(pause,mda+x-t1)
        }else pamt<-pause

       if(num.arm==1){
            sig21 <- sig31/rho21

            ### find design associated with rho  
            ### function supplied to root finder for step 5 
            ### check if uniroot has a solution (largest possible function value should be positive)
            tmp <- f51G2(t=t1+pamt,B.init=B.init,m.init=m.init,shape1=shape1,
                   scale1=scale1,x=x,n=n,mda=mda,l=l,pause=pause,
                   t1=t1,sig21=sig21
                   )
            if(tmp<0)  return(list(min=1e10))   
            t2.root <- try(uniroot(f51G2,interval=c(t1+pamt,mdaP+x),
                   B.init=B.init,m.init=m.init,shape1=shape1,
                   scale1=scale1,x=x,n=n,mda=mda,l=l,pause=pause,
                   t1=t1,sig21=sig21
                   ),silent=TRUE)
        }else{
            ## find design associated with rho11, rho21
            ### check if uniroot has a solution (largest possible function value should be positive)
            tmp <- f52G2(t=t1+pamt,B.init=B.init,m.init=m.init,shape1=shape1,
                   scale1=scale1,shape0=shape0,scale0=scale0,x=x,n=n,mda=mda,
                   l=l,r=r,pause=pause,t1=t1,sig31=sig31,sig30=sig30,
                   lam1=lam1,lam0=lam0,rho21=rho21
                   )
            if(tmp<0)  return(list(min=1e10))  
            t2.root <- try(uniroot(f52G2,interval=c(t1+pamt,mdaP+x),
                   B.init=B.init,m.init=m.init,shape1=shape1,
                   scale1=scale1,shape0=shape0,scale0=scale0,x=x,n=n,mda=mda,
                   l=l,r=r,pause=pause,t1=t1,sig31=sig31,sig30=sig30,
                   lam1=lam1,lam0=lam0,rho21=rho21
                   ),silent=TRUE)
        }

        if(class(t2.root)=="try-error")return(list(min=1e10)) 
        else if(suppressWarnings(warning()!="")){
          stop("the algorithm of uniroot() does not converge in 'maxiter' steps") 
        }else{ 
          t2 <- t2.root$root 
        }     

        ### mda adjusted for pause
        if(t2<mdaP)mdaP<-mdaP+pause 
        ### pause truncated for end of study
        if((t1>=mda) || (t2>=mda+pause)){
          pamt2<-min(pause,mda+x-t2)
       }else pamt2<-pause
      

        sig10 <- sqrt(sig2(B.init,m.init,shape0,scale0,x,n,t1+pamt,
                 mda,l,pause=pause,pt1=t1))
        sig20 <- sqrt(sig2(B.init,m.init,shape0,scale0,x,n,t2+pamt2,
                 mda,l,pause=pause,pt1=t1,pt2=t2))
        if(num.arm==2){
            sig11 <- sqrt(sig2(B.init,m.init,shape1,scale1,x,n,t1+pamt,
                  mda,l,pause=pause,pt1=t1))
            sig21 <- sqrt(sig2(B.init,m.init,shape1,scale1,x,n,t2+pamt2,
                  mda,l,pause=pause,pt1=t1,pt2=t2))
        }    

        rho10 <- sig30/sig10
        rho20 <- sig30/sig20
        if(rho10>=rho20-smallC) return(list(min=1e10))

        ## Step 7
        if(num.arm==1){
            u1 <- sqrt(n)*(log(lam0)-log(lam1))*lam1/sig11
            u2 <- sqrt(n)*(log(lam0)-log(lam1))*lam1/sig21
            u3 <- sqrt(n)*(log(lam0)-log(lam1))*lam1/sig31
            u <- c(u1,u2,u3)
        }else{
            v10 <- (sig10^2)/((1-r)*(lam0^2))
            v11 <- (sig11^2)/(r*(lam1^2))
            v20 <- (sig20^2)/((1-r)*(lam0^2))
            v21 <- (sig21^2)/(r*(lam1^2))
            v30 <- (sig30^2)/((1-r)*(lam0^2))
            v31 <- (sig31^2)/(r*(lam1^2))
            u1 <- sqrt(n)*(log(lam0)-log(lam1))/sqrt(v10+v11)
            u2 <- sqrt(n)*(log(lam0)-log(lam1))/sqrt(v20+v21) 
            u3 <- sqrt(n)*(log(lam0)-log(lam1))/sqrt(v30+v31)
            u <- c(u1,u2,u3)
        }
        # the number of patients accrued at the times of interim analyses
        # subtract pauses for compTime
        tinterim<-ifelse(t1>=mda,mda,t1)  ### remove pause from accrual total
        n1<-compTime(B.init,m.init,tinterim)
        ### remove pause from accrual total for unshifted compTime 
        tinterim<-ifelse(t2>=mdaP,mda,t2-pause)  
        n2<-compTime(B.init,m.init,tinterim)

        # the number of patients accrued at the times of interim analyses
        if(t1>=mda){
            n1<-compTime(B.init,m.init,mda)
            n2<-n1
        }else{
            n1<-compTime(B.init,m.init,t1)
            ### remove first pause from accrual total
            tinterim<-ifelse(t2>=mdaP,mda,t2-pause)  
            n2<-compTime(B.init,m.init,tinterim)
        }

        ## Steps 8-10
        sigma0 <- diag(3)
        sigma1 <- diag(3)
        sigma0[1,2] <- rho10/rho20
        sigma0[1,3] <- rho10
        sigma0[2,3] <- rho20
        sigma0[2,1] <- sigma0[1,2]
        sigma0[3,1] <- sigma0[1,3]
        sigma0[3,2] <- sigma0[2,3]
        
        sigma1[1,2] <- rho11/rho21
        sigma1[1,3] <- rho11
        sigma1[2,3] <- rho21
        sigma1[2,1] <- sigma1[1,2]
        sigma1[3,1] <- sigma1[1,3]
        sigma1[3,2] <- sigma1[2,3]
        
        alphaI1<-0
        alphaI2<-0
        if(sf=='futility'){
            C1U <- Inf   
            C2U <- Inf
            C3<-qnorm(1-alpha)
        }else{
            alphaI1<-spfun(rho10,alpha,sf)
            alphaI2<-spfun(rho20,alpha,sf)
            C1U <- qnorm(1-alphaI1)
            
            if(fC2U(0,C1U,sigma0,alphaI1,alphaI2)<0)return(list(min=1e10)) 
            C2Uroot <- try(uniroot(fC2U,c(0,qnorm(1-(alphaI2-alphaI1))),
                       C1U=C1U,sigma0=sigma0,alphaI1=alphaI1,alphaI2=alphaI2
            ),silent=TRUE)
            if(class(C2Uroot)=="try-error")return(list(min=1e10)) 
            if(suppressWarnings(warning()!="")){
               return(list(min=1e10))
            }else C2U <- C2Uroot$root

            if(fC3U(0,C1U,C2U,sigma0,alpha,alphaI2)<0)return(list(min=1e10)) 
            C3Uroot <- try(uniroot(fC3U,c(0,qnorm(1-(alpha-alphaI2))),
                           C1U=C1U,C2U=C2U,sigma0=sigma0,alpha=alpha,
                           alphaI2=alphaI2 
            ),silent=TRUE)
            if(class(C3Uroot)=="try-error"){
               return(list(min=1e10))
            }else C3 <- C3Uroot$root
        }


        if(C1L>C1U){
           return(list(min=1e10))  #invalid parameter region
        }else if(sign(powC2L(-3.5,C1L,C1U,C2U,C3,sigma1,beta,u1,u2,u3))==
                 sign(powC2L(qnorm(1-alpha),C1L,C1U,C2U,C3,sigma1,beta,u1,u2,u3))){
           return(list(min=1e10))  #when there is no solution
        }else{  ###C1 is C1L in the two-arm
           C2Lroot <- try(uniroot(powC2L,lower=-3.5,upper=qnorm(1-alpha),
                          C1L=C1L,C1U=C1U,C2U=C2U,C3=C3,sigma1=sigma1,
                          beta=beta,u1=u1,u2=u2,u3=u3
           ),silent=TRUE)  
           if(class(C2Lroot)=="try-error")return(list(min=1e10)) 
        }
        if(suppressWarnings(warning()!="")){
           return(list(min=1e10))  #when there is no solution
        }else if(C2Lroot$root>C2U){
           return(list(min=1e10))  #invalid parameter region
        }else{
           C2L <- C2Lroot$root
        }

        ps1Null<-pnorm(C1L)+(1-pnorm(C1U))
        ps3Null<-as.numeric(pmvnorm(lower=c(C1L,C2L),upper=c(C1U,C2U),
              sigma=sigma0[1:2,1:2]))
        ps2Null<-1-(ps1Null+ps3Null)
        ps1Alt<-pnorm(C1L-u1)+(1-pnorm(C1U-u1))
        ps3Alt<-as.numeric(pmvnorm(lower=c(C1L-u1,C2L-u2),upper=c(C1U-u1,C2U-u2),
              sigma=sigma0[1:2,1:2]))
        ps2Alt<-1-(ps1Alt+ps3Alt)

        EDA<- ifelse(t1<mdaP,ifelse(t2<mdaP,ps1Null*t1+
                ps2Null*t2+ps3Null*mdaP,
              ps1Null*t1+(1-ps1Null)*mdaP),mdaP)
        ETSL<-ps1Null*(t1+pamt)+
              ps2Null*(t2+pamt2)+
              ps3Null*(mdaP+x)
        ES<-  ps1Null*n1+ps2Null*n2+ps3Null*n
        pstopNull<-ps1Null+ps2Null
        pstopAlt<-ps1Alt+ps2Alt
        
    } ## end of three-stage algorithm
        
    
    if(target=="EDA")
      outcome1 <- EDA
    if(target=="ETSL")
      outcome1 <- ETSL
    if(target=="ES")
      outcome1<- ES
      
    if(num.stage==2){  
    return(list(min=outcome1,EDA=EDA,ETSL=ETSL,ES=ES,
            pstopNull=pstopNull,pstopAlt=pstopAlt,n1=n1,t1=t1,mda=mdaP,
            C1=C1,C1U=C1U,C2=C2,se=c(sig10,sig20,sig11,sig21),u=u))
    }
    if(num.stage==3){
    return(list(min=outcome1,EDA=EDA,ETSL=ETSL,ES=ES,
            pstopNull=pstopNull,pstopAlt=pstopAlt,n1=n1,n2=n2,t1=t1,
            t2=t2,mda=mdaP,
            C1L=C1L,C1U=C1U,C2L=C2L,C2U=C2U,C3=C3,se=c(sig10,sig20,
            sig30,sig11,sig21,sig31),u=u))
    }
  }

##############################################################################################



##############################################################################################

  ### wrapper function for f.Des so it returns the limited output 
  ### expected by function Optimize

  f.DesWrap <- function(optpar,B.init,m.init,alpha,beta,param,x,n,num.arm,r,num.stage,
                     target=c("EDA","ETSL","ES"),sf=c("futility","OF","Pocock"),
                     pause=0)
  {
       target <- match.arg(target)
       sf <- match.arg(sf)
       if(num.stage==3){
           rho1<-plogis(optpar[1])
           rho2<-plogis( optpar[1]+exp(optpar[2]) )
           optparin<-c(rho1,rho2,optpar[3])
       }else optparin<-optpar
       out<-f.Des(optpar=optparin,B.init=B.init,m.init=m.init,
              alpha=alpha,beta=beta,param=param,x=x,
              n=n,num.arm=num.arm,r=r,num.stage=num.stage,
              target=target,sf=sf,pause=pause) 
      return(out$min) 
  }
  
 
  

   

  ### numerical control values that can be changed by a user
  OptimDesControl<-function(trace=TRUE,tol=0.01,
                            n.int=c(1,5), 
                            aboveMin=c(1.05,1.10))
  {
  list(trace=trace,tol=tol,startC1L=-1.5,
       n.int=n.int,aboveMin=aboveMin)
  }

  ### alpha spending function
  spfun <- function(rho,alpha,sf=c("futility","OF","Pocock"))
  {
    sf <- match.arg(sf)
    if(sf=="futility")
      alpha.rho <- ifelse(rho<1&rho>=0,0,1)
    
    if(sf=="OF")
      alpha.rho <- min(2-2*pnorm(qnorm(1-alpha/2)/rho), alpha)
    
    if(sf=="Pocock")
      alpha.rho <- min(alpha*log(1+(exp(1)-1)*rho^2),alpha)
  
    return(alpha.rho)
  }

