np.OptimDes <-
function(B.init,m.init,alpha,beta,param,x,n=NULL,pn=NULL,pt=NULL,
         target=c("EDA","ETSL","ES"),sf=c("futility","OF","Pocock"),
         num.arm,r=0.5,num.stage=2,
         pause=0,
         control=OptimDesControl(), ...)
{
  ##set the seed for reproducibility of function optim
  if(exists('.Random.seed')){
    oldseed<-.Random.seed
  }else{
    set.seed(12357)
    oldseed<-.Random.seed
  }
  on.exit(set.seed(oldseed)) 

  
  target <- match.arg(target)
  sf <- match.arg(sf)

  if(!identical(num.arm,1) && !identical(num.arm,2))
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
  
  if(identical(num.arm,1)&!identical(sf,"futility"))
  {
    stop("only futility interim analysis for single-arm trial")
  }

  if(!num.stage%in%c(2,3))
  {
    stop("The optimal design allows two or three stages")
  }

  if(identical(num.arm,1)&&!identical(r,0.5))
  {
    warning("randomization ratio not equal to 0.5 for a single-arm study")
  }
  
  if(identical(num.arm,2)&& isTRUE(r*(1-r)<=0))
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

  if(!identical(target,"EDA") && !identical(target,"ETSL") && !identical(target,"ES"))
    stop("unknown target value")
  
  if(!identical(trace,FALSE) && !identical(trace,TRUE))
    stop("trace should be a logical value.")
        
  if(isTRUE(s(shape1,scale1,x)<=s(shape0,scale0,x)))
  {
    stop("the alternative survial rate should be greater than the null survival rate")
  }

  if(length(c(n,pn,pt))!=1)
    stop("There should be one input from n, pn or pt")

  b <- length(B.init)
  # single stage sample size, duration of accrual and study length
  fix.d <- FixDes(B.init,m.init,alpha,beta,param,x,num.arm,r)
  n0 <- fix.d$n0
  da <- fix.d$DA
  sl <- fix.d$SL

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

     nfix <- n0
     sl <- fix.d$SL

  if(!is.null(pn))
    n <- ceiling(pn*nfix)
  if(!is.null(pt))
  {
    t <- pt*sl-x
    n<-ceiling(compTime(B.init,m.init,t))
  }  
  
  if(identical(num.stage,2)){
        optout <- optimize(f=f.DesWrap,interval=c(0,1),tol=tol,B.init=B.init,
                m.init=m.init,alpha=alpha,beta=beta,param=param,x=x,n=n,
                num.arm=num.arm,r=r,num.stage=num.stage,target=target,sf=sf,
                pause=pause)          
       
        outcome <- optout$objective
        rho1 <- optout$minimum
        #outcome <- optout$objective
        #rho1 <- optout$minimum

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

       outcome <- optout$value
       rho1<-plogis(optout$par[1])
       rho2<-plogis( optout$par[1]+exp(optout$par[2]) )

       C1 <- optout$par[3]
       
       if(optout$convergence!=0 || optout$value==1e10){
         rho1 <- NA
         rho2 <- NA
         C1 <- NA 
       }     
    }
  
    if(outcome==1e10)stop('An optimal design could not be computed')
  
    ### re-run final converged values to obtain full output
    if(identical(num.stage,2)) {optpar <- rho1
    }else if(identical(num.stage,3)) optpar <-c(rho1,rho2,C1)  
  
    f.out <- f.Des(optpar,B.init,m.init,alpha,beta,param,x,n,
                   num.arm,r,num.stage,target,sf,pause=pause)
                   
  EDA <- f.out$EDA
  ETSL <- f.out$ETSL
  ES <- f.out$ES
  mda <- f.out$mda
  t1.last <- f.out$t1 
  n.last <- n
  n1 <- ceiling(f.out$n1)
  se<-f.out$se    
  u <-f.out$u
  pstopNull<-f.out$pstopNull
  pstopAlt<-f.out$pstopAlt
  C1U.last <- f.out$C1U
  if(identical(num.stage,2)){
    C1.last <- f.out$C1
    C2.last <- f.out$C2
  }else if(identical(num.stage,3)){
    t2.last<-f.out$t2
    n2 <- ceiling(f.out$n2)
    C1.last <- C1
    C2.last <- f.out$C2L
    C2U.last <- f.out$C2U
    C3.last <- f.out$C3
    }  
  
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

  if(identical(num.stage,2)){

   res <-structure(list(target=target,sf=sf,test=c(alpha=alpha,beta=beta,
            param=param,x=x), 
            design=c(num.arm=num.arm,r=r,num.stage=num.stage,pause=pause),
            accrual=list(B.init=B.init,m.init=m.init),
            result=c(EDA=EDA,ETSL=ETSL,ES=ES,pstopNull=pstopNull,pstopAlt=pstopAlt),
            n=c(n1=n1,n.last=n.last),
            stageTime=c(t1=t1.last,MTSL=mda+x,MDA=mda),boundary=c(C1L=C1.last,
            C1U=C1U.last,C2=C2.last),se=se,u=u,exposure=EW,all.info=NULL,
            single.stageTime=c(n0=n0,DA=da,SL=da+x,n0E=n0E,DAE=DAE,SLE=SLE)),
            class="OptimDes")
  }else if(identical(num.stage,3)){

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
            stageTime=c(t1=t1.last,t2=t2.last,MTSL=mda+x,MDA=mda),boundary=c(C1L=C1.last,
            C1U=C1U.last,C2L=C2.last,C2U=C2U.last,C3=C3.last),se=se,
            u=u,exposure=EW,all.info=NULL,
            single.stageTime=c(n0=n0,DA=da,SL=da+x,n0E=n0E,DAE=DAE,SLE=SLE)),
            class="OptimDes")
  }     
 
  return(res)
  
    
}


