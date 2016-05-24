"SimDes" <- 
function(object,B.init,m.init,weib0,weib1,interimRule='e1',
         sim.n=1000,e1conv=1/365,CMadj=F,attainI=1,attainT=1,FixDes="F",
         Rseed)
{  

  interimRule<-match.arg(interimRule,c('e1','n1','t1'))
  if(!identical(interimRule,'n1') && !identical(interimRule,'t1') && 
     !identical(interimRule,'e1'))stop('Invalid interimRule')


  if(missing(B.init))B.init<-object$accrual$B.init
  if(missing(m.init))m.init<-object$accrual$m.init
  if(!identical(length(B.init),length(m.init)))stop(
     "The length of the accrual rates and intervals must be equal")

  if(missing(weib0))weib0<-c(object$test["param1"],object$test["param2"])
  if(missing(weib1))weib1<-c(object$test["param3"],object$test["param4"])

  if(!missing(Rseed))set.seed(Rseed)

  if(any(identical(FixDes,"E"),identical(FixDes,"N")))interimRule<-'t1'

  cat("\n The interimRule is ",interimRule,"\n\n")

  num.arm<-as.numeric(object$design["num.arm"])
  num.stage<-as.numeric(object$design["num.stage"])
  x<-as.numeric(object$test["x"])
  alpha<-object$test["alpha"]
  pause<-as.numeric(object$design["pause"])
  C1L <- as.numeric(object$boundary[1])
  C1U <- as.numeric(object$boundary[2])
  if(identical(num.stage,2)){
    C2U <- as.numeric(object$boundary[3])
  }else{
    C2L<-as.numeric(object$boundary[3])
    C2U<-as.numeric(object$boundary[4])
    C3U<-as.numeric(object$boundary[5])
  }
  r<-as.numeric(object$design["r"])
  futility<-identical(as.character(object$sf),'futility')
  if(CMadj && !identical(num.arm,1))stop("CM adjustment available with num.arm=1 only")
  if(isTRUE(CMadj) && identical(num.stage,3))  
    stop("The CM adjustment is for two-stage trials only.")  


  if(!identical(length(B.init),length(m.init)))
  {
    stop("projected patient times and numbers should be of equal length")
  }

  if(isTRUE(x>=max(B.init)))
  {
    stop("the survival time of interest is beyond the projected accrual time")
  } 

  ### apply CM correction factor

  Stime<-object$single.stageTime

  ### apply adjustment to exact fixed design size
  ### per CM 2003
  if (CMadj){
        CMadjfac=as.numeric(Stime["n0E"]/Stime["n0"])
        exposure<-as.numeric(object$exposure[2])
  }else{
        CMadjfac<-1
        exposure<-as.numeric(object$exposure[1])
        if(identical(num.stage,3)){
            exposure2<-as.numeric(object$exposure[2])
        }
  }

  if(identical(FixDes,"E")){n<-as.numeric(Stime["n0E"])
  }else if(identical(FixDes,"N")){
    n<-as.numeric(Stime["n0"])
    cnorm<-qnorm(1-alpha)
  }else{
      ### the attained interim sample sizes are re-computed at 
      ### the beginning of each
      ### iteration because they may change based on interimRule.
      ### they are computed here to check consistency with overall
      ### sample size
      if(identical(num.stage,2)){
          n<-ceiling(CMadjfac*attainT*as.numeric(object$n[2]))
          if(isTRUE(n>sum(m.init)))stop('The attained sample size exceeds the maximum accrual')
          n1<-ceiling(CMadjfac*attainI*as.numeric(object$n[1]))
          if(isTRUE(n1>n))stop("Single-stage design is optimal")
      }else{
          n<-ceiling(attainT*as.numeric(object$n[3]))
          if(isTRUE(n>sum(m.init)))stop('Projected interim size exceeds final size')
          n1<-ceiling(attainI*as.numeric(object$n[1]))
          n2<-ceiling(attainI*as.numeric(object$n[2]))
          if(isTRUE(n2>n))stop("Projected interim size exceeds final size")
      }
  }
  ### if interim analyses are planned after the end of accrual
  ### do not allow option interimRule='n1'
  if((identical(interimRule,'n1')  && isTRUE(n1>=n)) || (identical(num.stage,3) && n2>=n))stop(
     'interimRule n1 not allowed when the interim analysis is scheduled after the end of accrual')

  if(identical(num.arm,2))ntrt<-ceiling(r*n)
  else ntrt<- n
  ncntl<- n-ntrt

  shape0 <- weib0[1]
  scale0 <- weib0[2]
  shape1 <- weib1[1]
  scale1 <- weib1[2]  

  ### event rate under null hypothesis (not event-free as in OptimDes)
  p0<-pweibull(x,shape=shape0,scale=scale0)
  
  ### cutoff for number of events if exact binomial test is used at second stage for one-arm design
  if(identical(num.arm,1)) {
      cexact<-qbinom(alpha,n,p0)       # one-arm
      if(pbinom(cexact,n,p0)>alpha)cexact<-cexact-1
  }

  ### compute mda and distribution restricted to (0,mda)
  ### the CMadjfac is applied to mda through the adjustment to n
  ### times are adjusted for accrual pauses later
  mdaout<-compMDA(B.init,m.init,n)
  mdaU<-mdaout[1]
  l<-mdaout[2]

  B <- B.init[1:l]  # "l" not "one"
  B[l]<-mdaU
  M <- m.init[1:l]
  M[l]<- ifelse(l>1, n-sum(m.init[1:(l-1)]) , n )
  Pi<-M/sum(M)  # proportion in each accrual interval
  B<-c(0,B)     # add B0 for diff calculation

  ### bounds needed to compute truncated exposure conditional on Y<=t1
  tlow<-0
  thigh<-max(B.init)+x+2*pause
  
  alphaNorm <- 0
  alphaExact<- 0
  powerNorm   <- 0
  powerExact  <- 0
  eda<-0
  etsl<-0
  es<-0
  edaAlt<-0
  etslAlt<-0
  esAlt<-0
  aveE<-0
  aveE2<-0
  intvarNull<-0
  intvarNull2<-0
  intvarNull3<-0
  intvarAlt<-0
  intvarAlt2<-0
  intvarAlt3<-0
  n1int  <- 0
  t1int <- 0
  n2int  <- 0
  t2int <- 0
  difIntFutL<- Inf
  difIntFutH<- -Inf
  difFinSupL<- Inf
  difFinFutH<- -Inf
  difIntSupL<- Inf
  difIntSupH<- -Inf
  lambda0<-lambda(shape0,scale0,x)
  pstopNull<-0
  pstopAlt<-0
  pstopENull<-0
  pstopEAlt<-0

  for(i in 1:sim.n)
  {
   
    ### create study entry time for each patient without accrual pauses
    randI <- sample(1:l,ntrt,replace=TRUE,prob=Pi) #sample accrual intervals
    u<-runif(ntrt)
    sample.Y1 <- sort( B[randI] + u*(B[randI+1]-B[randI]) )
    if(identical(num.arm,1)) {
        sample.Y0<-NULL
    }else{
        randI <- sample(1:l,ncntl,replace=TRUE,prob=Pi) #sample accrual intervals
        u<-runif(ncntl)
        sample.Y0 <- sort( B[randI] + u*(B[randI+1]-B[randI]) )
    }
    sample.Y<-c(sample.Y1,sample.Y0)
    ### reset mda unadjusted for pauses to value obtained
    mdaU<-max(sample.Y)
   

    # for 'e1' and 't1', the n are derived later
    ### compute time for interim to match planned exposure
    if(identical(interimRule,'e1')){  
       ### exposure and mda already have the CMadj adjustment 
       t1.root <- try(uniroot(ft,c(tlow,thigh),tol=e1conv,y=sample.Y,
                  x=x,EW=exposure,pause=pause,mda=mdaU))

       if(class(t1.root)=="try-error"){
            stop("Function uniroot did not converge for exposure time calculation")
       }else if(suppressWarnings(!identical(warning(),"")))
          stop("Function uniroot did not converge for exposure time calculation")
       else{
          t1<-t1.root$root
       }
       t1<-t1*attainI
       if(t1<x+sample.Y1[1]-pause){
          warning("Matched exposure yielded interim time less than x from enrollment")
          t1<-x+sample.Y1[1]-pause
       }
       if(identical(num.arm,2)){
           if(t1<x+sample.Y0[1]-pause){
              warning("Matched exposure yielded interim time less than x from enrollment")
              t1<-x+sample.Y0[1]-pause
           }
       }
       if(identical(num.stage,3)){
           t2.root <- try(uniroot(ft,c(tlow,thigh),tol=e1conv,y=sample.Y,
                      x=x,EW=exposure2,pause=pause,pt1=t1,mda=mdaU))
           if(class(t1.root)=="try-error"){
                stop("Function uniroot did not converge for exposure time calculation")
           }else if(suppressWarnings(!identical(warning(),"")))
                stop("Function uniroot did not converge for exposure time calculation")
           else{
               t2<-t2.root$root
           }
           t2<-t2*attainI
           if(t2<x+sample.Y1[1]-pause){
              warning("Matched exposure yielded second interim time less than x from enrollment")
              t2<-x+sample.Y1[1]-pause
           }
           if(identical(num.arm,2)){
               if(t2<x+sample.Y0[1]-pause){
                  warning("Matched exposure yielded second interim time less than x from enrollment")
                  t2<-x+sample.Y0[1]-pause
               }
           }
        }
    }else if(identical(interimRule,'t1')){
       t1<-object$stageTime[1]*CMadjfac*attainI
       if(isTRUE(t1<x+sample.Y1[1]-pause)){
          warning("Adjusted t1 occurred before time x from enrollment")
          t1<-x+sample.Y1[1]-pause
       }
       if(identical(num.arm,2)){
           if(isTRUE(t1<x+sample.Y0[1]-pause)){
              warning("Adjusted t1 occurred before time x from enrollment")
              t1<-max(t1,x+sample.Y0[1]-pause)
           }
       }
       if(identical(num.stage,3)){
           t2<-object$stageTime[2]**attainI
           if(isTRUE(t2<x+sample.Y1[1]-pause)){
              warning("Adjusted t2 occurred before time x from enrollment")
              t2<-x+sample.Y1[1]-pause
           }
           if(identical(num.arm,2)){
               if(isTRUE(t2<x+sample.Y0[1]-pause)){
                  warning("Adjusted t2 occurred before time x from enrollment")
                  t2<-max(t2,x+sample.Y0[1]-pause)
               }
           }
       }
    }else{
       # the interim times are set later when rule 'n1' is used
       n1<-ceiling(as.numeric(object$n[1])*CMadjfac*attainI)
       if(identical(num.stage,3))n2<-ceiling(as.numeric(object$n[2])*attainI)
       ntrt1<-n1
       ncntl1<-0
       if(identical(num.arm,2)){
          ntrt1<-rbinom(1,n1,ntrt/n)
          ncntl1<-n1-ntrt1
       }
       if(identical(num.stage,3)){
           ntrt2<-n2
           ncntl2<-0
           if(identical(num.arm,2)){
              ntrt2<-ntrt1+rbinom(1,n2-n1,ntrt/n)
              ncntl2<-n2-ntrt2
           }
       }
    }  

    ### deterimine mda with pauses
    if(identical(interimRule,'n1')){
        ### both interims must be before mda
        if(identical(num.stage,2))mda<-mdaU+pause else mda<-mdaU+2*pause
    }else{
        if(t1<mdaU)mda<-mdaU+pause else mda<-mdaU
        if(identical(num.stage,3) && (t2<mda))mda<-mda+pause
    }

    ### check potential inconsistency due to high attained time for interim
    if(!identical(interimRule,'n1')){
        if(isTRUE(t1>mda+x))stop('The targetted interim analysis was scheduled after the final analysis')
        if(isTRUE(identical(num.stage,3))){ 
            if(t2>mda+x)stop('The targetted interim analysis was scheduled after the final analysis')
        }
    }

    ### add pause(s) to accrual times
    if(!identical(interimRule,'n1')){
        if(t1<mda)pamt<-pause else pamt<-0
        sample.Y[sample.Y>t1]<-sample.Y[sample.Y>t1]+pamt
        sample.Y1[sample.Y1>t1]<-sample.Y1[sample.Y1>t1]+pamt
        if(identical(num.arm,2))sample.Y0[sample.Y0>t1]<-sample.Y0[sample.Y0>t1]+pamt
        if(identical(num.stage,3)){
            if(t2<mda)pamt2<-pause else pamt2<-0
            sample.Y[sample.Y>t2]<-sample.Y[sample.Y>t2]+pamt2
            sample.Y1[sample.Y1>t2]<-sample.Y1[sample.Y1>t2]+pamt2
            if(identical(num.arm,2))sample.Y0[sample.Y0>t2]<-sample.Y0[sample.Y0>t2]+pamt2
        }
    }else{  ## interims based on sample sizes, interims cannot occur after mda
        sample.Y[(n1+1):n]<-sample.Y[(n1+1):n]+pause
        sample.Y1[(ntrt1+1):ntrt]<-sample.Y1[(ntrt1+1):ntrt]+pause
        if(identical(num.arm,2)){
            sample.Y0[(ncntl1+1):ncntl]<-sample.Y0[(ncntl1+1):ncntl]+pause
        }
        if(identical(num.stage,3)){
            sample.Y[(n2+1):n]<-sample.Y[(n2+1):n]+pause
            sample.Y1[(ntrt2+1):ntrt]<-sample.Y1[(ntrt2+1):ntrt]+pause
            if(identical(num.arm,2)){
                sample.Y0[(ncntl2+1):ncntl]<-sample.Y0[(ncntl2+1):ncntl]+pause
            }
        }
    }


    ### alpha level

    ### event times
    sample.T1 <- rweibull(ntrt,shape=shape0,scale=scale0) #null survival times
    if(identical(num.arm,2)){
        sample.T0 <- rweibull(ncntl,shape=shape0,scale=scale0) 
    }else sample.T0<-NULL

    ### exact analyses based on full data set
    obs1<-sum(sample.T1<=x)
    if(identical(num.arm,2)){
       obs0<-sum(sample.T0<=x)
       pv <- fisher.test(matrix(c(obs1,obs0,ntrt-obs1,ncntl-obs0),nrow=2),
               alternative="less",conf.int=FALSE)$p.value
    }

    if(identical(FixDes,"E")){
        if(identical(num.arm,1)){
            alphaExact<- alphaExact + (obs1<=cexact)
        }else alphaExact<- alphaExact + (pv<=alpha)
    }else if(identical(FixDes,"N")){
         outN <- TestStage(mda+x,1,x, num.arm, num.stage, 
                   sample.Y1, sample.T1, sample.Y0, sample.T0, 
                   p0, C1L, C1U,    ## only Z-stat is used
                   printTest=FALSE)   
         z<-as.numeric(outN["z"])
         if(z>=cnorm)alphaNorm<-alphaNorm+1
    }else{  ### multi-stage designs

        ### create accrual and event time variables for each analysis stage
        if(any(identical(interimRule,'t1'),identical(interimRule,'e1'))){
           ntrt1<-sum(sample.Y1<=t1)
           ncntl1<-0
           if(identical(num.arm,2))ncntl1 <- sum(sample.Y0<=t1)
           n1<-ntrt1+ncntl1
           if(identical(num.stage,3)){
               ntrt2<-ntrt1+sum(sample.Y1>t1 & sample.Y1<=t2)
               ncntl2<-0
               if(identical(num.arm,2))ncntl2 <- ncntl1+sum(sample.Y0>t1 & sample.Y0<=t2)
               n2<-ntrt2+ncntl2
           }
           Y11 <- sample.Y1[sample.Y1<=t1]
           T11 <- sample.T1[sample.Y1<=t1]
           if(identical(num.stage,2)){
               Y21 <- sample.Y1[sample.Y1>t1]
               T21 <- sample.T1[sample.Y1>t1]
           }else{
               Y21 <- sample.Y1[sample.Y1>t1 & sample.Y1<=t2]
               T21 <- sample.T1[sample.Y1>t1 & sample.Y1<=t2]
               Y31 <- sample.Y1[sample.Y1>t2]
               T31 <- sample.T1[sample.Y1>t2]
           }
           if(identical(num.arm,2)) {
               Y10 <- sample.Y0[sample.Y0<=t1]
               T10 <- sample.T0[sample.Y0<=t1]
               if(identical(num.stage,2)){
                   Y20 <- sample.Y0[sample.Y0>t1]
                   T20 <- sample.T0[sample.Y0>t1]
               }else{
                   Y20 <- sample.Y0[sample.Y0>t1 & sample.Y0<=t2]
                   T20 <- sample.T0[sample.Y0>t1 & sample.Y0<=t2]
                   Y30 <- sample.Y0[sample.Y0>t2]
                   T30 <- sample.T0[sample.Y0>t2]
               }
           }else{
               Y10 <- NULL
               T10 <- NULL
               Y20 <- NULL
               T20 <- NULL
               Y30 <- NULL
               T30 <- NULL
           }
        }else{       ### interim based on target sample size unless
                     ### all patients are accrued before the interim
           if(isTRUE(x>pause+sample.Y1[ntrt1]-sample.Y1[1])){  ## interim too early
              t1<-x+sample.Y1[1]-pause
              ntrt1<-min(ntrt,sum(sample.Y1<t1)+1)
              t1<-sample.Y1[ntrt1]
              if(identical(num.stage,3)){
                ntrt2<- max(ntrt1,ntrt2)
                t2<-max(sample.Y1[ntrt2],t1)
              }
              warning("ntrt1 patients were accrued before time x from enrollment")
           }else{
              t1<-sample.Y1[ntrt1]
              if(identical(num.stage,3))t2<-sample.Y1[ntrt2]
           }

           if(identical(num.arm,2)){  ## interim too early
               if(isTRUE(x>pause+sample.Y0[ncntl1]-sample.Y0[1])){
                   t1<-max(t1,x+sample.Y0[1]-pause)
                   ncntl1<-min(ncntl,sum(sample.Y0<t1)+1)
                   t1<-max(t1,sample.Y0[ncntl1])
                   if(identical(num.stage,3)){
                      ncntl2<- max(ncntl1,ncntl2)
                      t2<-max(t1,t2,sample.Y0[ncntl2])
                   }
                   warning("ncntl1 patients were accrued before time x from enrollment")
               }else{  
                  t1<-max(t1,sample.Y0[ncntl1])
                  if(identical(num.stage,3))t2<-max(t2,sample.Y0[ncntl2])
               }
           }

           Y11 <- sample.Y1[1:ntrt1]
           T11 <- sample.T1[1:ntrt1]
           if(identical(num.stage,2)){
               if(isTRUE(ntrt1<ntrt)){
                  Y21 <- sample.Y1[(ntrt1+1):ntrt]
                  T21 <- sample.T1[(ntrt1+1):ntrt]
               }else{
                  Y21<-NULL
                  T21<-NULL
               }
           }else{
               Y21<-NULL
               T21<-NULL
               Y31<-NULL
               T31<-NULL
               if(isTRUE(ntrt1<ntrt2)){
                  Y21 <- sample.Y1[(ntrt1+1):ntrt2]
                  T21 <- sample.T1[(ntrt1+1):ntrt2]
               }
               if(isTRUE(ntrt2<ntrt)){
                  Y31 <- sample.Y1[(ntrt2+1):ntrt]
                  T31 <- sample.T1[(ntrt2+1):ntrt]
               }
           }
           
           if(identical(num.arm,2)) {
               Y10 <- sample.Y0[1:ncntl1]
               T10 <- sample.T0[1:ncntl1]
               if(identical(num.stage,2)){
                   if(isTRUE(ncntl1<ncntl)){
                      Y20 <- sample.Y0[(ncntl1+1):ncntl]
                      T20 <- sample.T0[(ncntl1+1):ncntl]
                   }else{
                      Y20<-NULL
                      T20<-NULL
                   }
               }else{
                   Y20<-NULL
                   T20<-NULL
                   Y30<-NULL
                   T30<-NULL
                   if(isTRUE(ncntl1<ncntl2)){
                      Y20 <- sample.Y0[(ncntl1+1):ncntl2]
                      T20 <- sample.T0[(ncntl1+1):ncntl2]
                   }
                   if(isTRUE(ncntl2<ncntl)){
                      Y30 <- sample.Y0[(ncntl2+1):ncntl]
                      T30 <- sample.T0[(ncntl2+1):ncntl]
                   }
               }
           }else{
               Y10 <- NULL
               T10 <- NULL
               Y20 <- NULL
               T20 <- NULL
               Y30 <- NULL
               T30 <- NULL
           }
        }


        ### set pauses between t1, t2 and times of analyses
        ### note that pamt is pause until analysis here and not
        ### pause until accrual resumes, which is never if 
        ### the analysis occurs after mda
        if(t1<mda)pamt<-pause else pamt<-min(pause,max(mda+x-t1,0))
        if(identical(num.stage,2)){
            pamt2<-0 
        }else{
           if(t2<mda)pamt2<-pause else pamt2<-min(pause,max(mda+x-t2,0)) 
        }
        
        ### compute interim and z-statistics
        ### the final z-statistic is computed even if the interim
        ### is less than c1.  the stopping rule is imposed below

         out1 <- TestStage(t1+pamt,1,x, num.arm, num.stage, 
                   Y11, T11, Y10, T10, p0, C1L, C1U,    
                   printTest=FALSE)
         
         T1<-c(T11,T21)
         T0<-c(T10,T20)
         Y1<-c(Y11,Y21)
         Y0<-c(Y10,Y20)

         if(identical(num.stage,3)){
           T1<-c(T1,T31)
           T0<-c(T0,T30)
           Y1<-c(Y1,Y31)
           Y0<-c(Y0,Y30)
         }

         tint2<-ifelse(identical(num.stage,2),mda+x,t2+pamt2)
         out2 <- TestStage(tint2,2,x, num.arm, num.stage, 
                   Y1, T1, Y0, T0, p0, C2L=C2L, C2U=C2U,    
                   printTest=FALSE)

         if(identical(num.stage,3)){
             out3 <- TestStage(mda+x,3,x, num.arm, num.stage, 
                       Y1, T1, Y0, T0, p0, C3U=C3U, 
                       printTest=FALSE)
         }

        Zint1<-as.numeric(out1["z"])
        Zint2<-as.numeric(out2["z"])
        if(identical(num.stage,3))Zint3<-as.numeric(out3["z"])
        

        int1<-isTRUE(Zint1<=C1L || Zint1>=C1U)
        ### accumulate planned exposure and second interim time
        ### even if trial stopped at first interim
        n1int<-n1int+n1
        t1int<-t1int+t1
        aveE<-aveE + sum(pmin(x,pmax(0,t1+pamt-sample.Y)))
        if((identical(num.stage,3))){
            aveE2<-aveE2 + 
                   sum(pmin(x,pmax(0,t2+pamt2-sample.Y)))
            n2int<-n2int+n2
            t2int<-t2int+t2
            t2m<-t2+pamt2
            n2m<-n2
        }

        if(identical(num.stage,3) && !int1){
            int2<-isTRUE(Zint2<=C2L || Zint2>=C2U)
        }else {
            int2<-0
            t2m<-0      ### arbitrary t2 and n2 values for 2-stage designs
            n2m<-0      ### to simplify calculations
        }
        eda<-eda + (!int1 && !int2)*max(sample.Y1,sample.Y0) + 
                   int1*max(Y11,Y10) + int2*max(Y11,Y10,Y21,Y20)
        etsl<-etsl + (!int1 && !int2)*(max(sample.Y1,sample.Y0)+x) + 
                      int1*(t1+pamt) + int2*t2m
        es<-es + (!int1 && !int2)*n + int1*n1 + int2*n2m
        pstopNull<-pstopNull + int1*(n1<n) + int2*(n2m<n)
  
        if(!int1 && !int2) {
           if(identical(num.arm,1))alphaExact<- alphaExact + (obs1<=cexact)
           else alphaExact<- alphaExact + (pv<=alpha)
        }


        dif1<-exp(-exp(out1["Log(cumL1)"])) - exp(-exp(out1["Log(cumL0)"]))
        dif2<-exp(-exp(out2["Log(cumL1)"])) - exp(-exp(out2["Log(cumL0)"]))
        if(identical(num.stage,3))dif3<-exp(-exp(out3["Log(cumL1)"])) - exp(-exp(out3["Log(cumL0)"]))
  
        ##  l low,  h high
        ##  R final planned analysis
        ##  I interim superiority
        ##  K interim futility
        if(identical(num.stage,2)){
            # stop for efficacy at interim
            if(isTRUE(Zint1>=C1U)){
                alphaNorm<-alphaNorm+1
                difIntSupL<-min(difIntSupL,dif1)
                pstopENull<-pstopENull+1
            # effiacy at final analysis
            }else if(isTRUE(Zint1>C1L)&& isTRUE(Zint2>=C2U)){
              alphaNorm <- alphaNorm + 1
              difFinSupL<-min(difFinSupL,dif2)
              difIntFutL<-min(difIntFutL,dif1)
              difIntSupH<-max(difIntSupH,dif1)
            # futility at final
            }else if(isTRUE(Zint1>C1L)){
              difFinFutH<-max(difFinFutH,dif2)
              difIntFutL<-min(difIntFutL,dif1)
              difIntSupH<-max(difIntSupH,dif1)
            # futility at interim
            }else{
               difIntFutH<-max(difIntFutH,dif1)
               difIntSupH<-max(difIntSupH,dif1)
            }
        }else{ # num.stage is 3
            # efficacy at first interim
            if(isTRUE(Zint1>=C1U)){
                alphaNorm<-alphaNorm+1
                difIntSupL<-min(difIntSupL,dif1)
                pstopENull<-pstopENull+1
            # efficacy at the second interim
            }else if(isTRUE(Zint1>C1L)&& isTRUE(Zint2>=C2U)){
              alphaNorm <- alphaNorm + 1
              difIntSupL<-min(difIntSupL,dif2)
              difIntSupH<-max(difIntSupH,dif1)
              pstopENull<-pstopENull+1
            # effiacy at final
            }else if(isTRUE(Zint1>C1L)&& isTRUE(Zint2>C2L) && isTRUE(Zint3>=C3U)){
              alphaNorm <- alphaNorm + 1
              difFinSupL<-min(difFinSupL,dif3)
              difIntFutL<-min(difIntFutL,dif1,dif2)
              difIntSupH<-max(difIntSupH,dif1,dif2)
            # futlity at final
            }else if(isTRUE(Zint1>C1L) && isTRUE(Zint2>C2L)){
              difIntFutL<-min(difIntFutL,dif1,dif2)
              difIntSupH<-max(difIntSupH,dif1,dif2)
              difFinFutH<-max(difFinFutH,dif3)
            # futlity at second interim
            }else if(isTRUE(Zint1>C1L)){
              difIntSupH<-max(difIntSupH,dif1)
              difIntFutL<-min(difIntFutL,dif1)
              difIntFutH<-max(difIntFutH,dif2)
            # futility at first interim
            }else{
               difIntFutH<-max(difIntFutH,dif1)
            }
        }

        ### proportion of info at interim
        intvarNull<-intvarNull+(out1["se"])^2
        intvarNull2<-intvarNull2+(out2["se"])^2
        if(identical(num.stage,3))intvarNull3<-intvarNull3+(out3["se"])^2
    } 
     
     

    ### power
    ### re-use patient start times and associated calculations
    sample.T1 <- rweibull(ntrt,shape=shape1,scale=scale1) #alt survival times

    ### exact analyses based on full data set
    obs1<-sum(sample.T1<=x)
    if(identical(num.arm,2)){
       obs0<-sum(sample.T0<=x)
       pv <- fisher.test(matrix(c(obs1,obs0,ntrt-obs1,ncntl-obs0),nrow=2),
               alternative="less",conf.int=FALSE)$p.value
    }

    if(identical(FixDes,"E")){
        if(identical(num.arm,1)){
            powerExact<- powerExact + (obs1<=cexact)
        }else powerExact<- powerExact + (pv<=alpha)
    }else if(identical(FixDes,"N")){
         outN <- TestStage(mda+x,1,x, num.arm, num.stage, 
                   sample.Y1, sample.T1, sample.Y0, sample.T0, 
                   p0, C1L, C1U,    ## only Z-stat is used
                   printTest=FALSE)   
         z<-as.numeric(outN["z"])
         if(z>=cnorm)powerNorm<-powerNorm+1
    }else{  ### multi-stage designs

        ### create accrual and event time variables for each analysis stage
        if(any(identical(interimRule,'t1'),identical(interimRule,'e1'))){
           T11 <- sample.T1[sample.Y1<=t1]
           if(identical(num.stage,2)){
               T21 <- sample.T1[sample.Y1>t1]
           }else{
               T21 <- sample.T1[sample.Y1>t1 & sample.Y1<=t2]
               T31 <- sample.T1[sample.Y1>t2]
           }
        }else{       ### interim based on target sample size unless
                     ### all patients are accrued before the interim

           T11 <- sample.T1[1:ntrt1]
           if(identical(num.stage,2)){
               if(isTRUE(ntrt1<ntrt)){
                  T21 <- sample.T1[(ntrt1+1):ntrt]
               }else{
                  T21<-NULL
               }
           }else{
               T21<-NULL
               T31<-NULL
               if(isTRUE(ntrt1<ntrt2)){
                  T21 <- sample.T1[(ntrt1+1):ntrt2]
               }
               if(isTRUE(ntrt2<ntrt)){
                  T31 <- sample.T1[(ntrt2+1):ntrt]
               }
           }
        }

        
        ### compute interim and z-statistics
        ### the final z-statistic is computed even if the interim
        ### is less than c1.  the stopping rule is imposed below

         out1 <- TestStage(t1+pamt,1,x, num.arm, num.stage, 
                   Y11, T11, Y10, T10, p0, C1L, C1U,    
                   printTest=FALSE)
         
         T1<-c(T11,T21)
         if(identical(num.stage,3)){
           T1<-c(T1,T31)
         }

         tint2<-ifelse(identical(num.stage,2),mda+x,t2+pamt2)
         out2 <- TestStage(tint2,2,x, num.arm, num.stage, 
                   Y1, T1, Y0, T0, p0, C2L=C2L, C2U=C2U,    
                   printTest=FALSE)

         if(identical(num.stage,3)){
             out3 <- TestStage(mda+x,3,x, num.arm, num.stage, 
                       Y1, T1, Y0, T0, p0, C3U=C3U, 
                       printTest=FALSE)
         }

        Zint1<-as.numeric(out1["z"])
        Zint2<-as.numeric(out2["z"])
        if(identical(num.stage,3))Zint3<-as.numeric(out3["z"])
        
        int1<-isTRUE(Zint1<=C1L || Zint1>=C1U)

        if(identical(num.stage,3)){
            t2m<-t2+pamt2
            n2m<-n2
        }

        if(identical(num.stage,3) && !int1){
            int2<-isTRUE(Zint2<=C2L || Zint2>=C2U)
        }else {
            int2<-0
            t2m<-0      ### arbitrary t2 and n2 values for 2-stage designs
            n2m<-0      ### to simplify calculations
        }
        edaAlt<-edaAlt + (!int1 && !int2)*max(sample.Y1,sample.Y0) + 
                   int1*max(Y11,Y10) + int2*max(Y11,Y10,Y21,Y20)
        etslAlt<-etslAlt + (!int1 && !int2)*(max(sample.Y1,sample.Y0)+x) + 
                      int1*t1 + int2*t2m
        esAlt<-esAlt + (!int1 && !int2)*n + int1*n1 + int2*n2m

        pstopAlt<-pstopAlt + int1*(n1<n) + int2*(n2m<n)

        if(!int1 && !int2) {
           if(identical(num.arm,1))powerExact<- powerExact + (obs1<=cexact)
           else powerExact<- powerExact + (pv<=object$test["alpha"])
        }


        dif1<-exp(-exp(out1["Log(cumL1)"])) - exp(-exp(out1["Log(cumL0)"]))
        dif2<-exp(-exp(out2["Log(cumL1)"])) - exp(-exp(out2["Log(cumL0)"]))
        if(identical(num.stage,3))dif3<-exp(-exp(out3["Log(cumL1)"])) - exp(-exp(out3["Log(cumL0)"]))

        ##  l low,  h high
        ##  R final planned analysis
        ##  I interim superiority
        ##  K interim futility
        if(identical(num.stage,2)){
            # stop for efficacy at interim
            if(isTRUE(Zint1>=C1U)){
                powerNorm<-powerNorm+1
                difIntSupL<-min(difIntSupL,dif1)
                pstopEAlt<-pstopEAlt+1
            # effiacy at final analysis
            }else if(isTRUE(Zint1>C1L)&& isTRUE(Zint2>=C2U)){
              powerNorm <- powerNorm + 1
              difFinSupL<-min(difFinSupL,dif2)
              difIntFutL<-min(difIntFutL,dif1)
              difIntSupH<-max(difIntSupH,dif1)
            # futility at final
            }else if(isTRUE(Zint1>C1L)){
              difFinFutH<-max(difFinFutH,dif2)
              difIntFutL<-min(difIntFutL,dif1)
              difIntSupH<-max(difIntSupH,dif1)
            # futility at interim
            }else{
               difIntFutH<-max(difIntFutH,dif1)
               difIntSupH<-max(difIntSupH,dif1)
            }
        }else{ # num.stage is 3
            # efficacy at first interim
            if(isTRUE(Zint1>=C1U)){
                powerNorm<-powerNorm+1
                difIntSupL<-min(difIntSupL,dif1)
                pstopEAlt<-pstopEAlt+1
            # efficacy at the second interim
            }else if(isTRUE(Zint1>C1L)&& isTRUE(Zint2>=C2U)){
              powerNorm <- powerNorm + 1
              difIntSupL<-min(difIntSupL,dif2)
              difIntSupH<-max(difIntSupH,dif1)
              pstopEAlt<-pstopEAlt+1
            # effiacy at final
            }else if(isTRUE(Zint1>C1L)&& isTRUE(Zint2>C2L) && isTRUE(Zint3>=C3U)){
              powerNorm <- powerNorm + 1
              difFinSupL<-min(difFinSupL,dif3)
              difIntFutL<-min(difIntFutL,dif1,dif2)
              difIntSupH<-max(difIntSupH,dif1,dif2)
            # futlity at final
            }else if(isTRUE(Zint1>C1L) && isTRUE(Zint2>C2L)){
              difIntFutL<-min(difIntFutL,dif1,dif2)
              difIntSupH<-max(difIntSupH,dif1,dif2)
              difFinFutH<-max(difFinFutH,dif3)
            # futlity at second interim
            }else if(isTRUE(Zint1>C1L)){
              difIntSupH<-max(difIntSupH,dif1)
              difIntFutL<-min(difIntFutL,dif1)
              difIntFutH<-max(difIntFutH,dif2)
            # futility at first interim
            }else{
               difIntFutH<-max(difIntFutH,dif1)
            }
        }

        ### proportion of info at interim
        intvarAlt<-intvarAlt+(out1["se"])^2
        intvarAlt2<-intvarAlt2+(out2["se"])^2
        if(identical(num.stage,3))intvarAlt3<-intvarAlt3+(out3["se"])^2
    }
  }


  alphaExact <- alphaExact/sim.n
  alphaNorm  <- alphaNorm/sim.n
  if(!futility){alphaExact<-NA; powerExact<-NA}
  powerExact   <- powerExact/sim.n
  powerNorm    <- powerNorm/sim.n
  eda   <- eda/sim.n
  etsl  <- etsl/sim.n
  es    <- es/sim.n
  edaAlt   <- edaAlt/sim.n
  etslAlt  <- etslAlt/sim.n
  esAlt    <- esAlt/sim.n
  aveE  <- aveE/sim.n
  aveE2 <- aveE2/sim.n
  n1int <- n1int/sim.n
  t1int <- t1int/sim.n
  n2int <- n2int/sim.n
  t2int <- t2int/sim.n
  intvarNull<-intvarNull/sim.n
  intvarNull2<-intvarNull2/sim.n
  intvarNull3<-intvarNull3/sim.n
  finvarNull<-ifelse(identical(num.stage,2),intvarNull2,intvarNull3)
  pinfoNull <- finvarNull/intvarNull
  if(identical(num.stage,3))pinfoNull2 <- finvarNull/intvarNull2
  pstopNull <- pstopNull/sim.n
  pstopENull <- pstopENull/sim.n
  intvarAlt<-intvarAlt/sim.n
  intvarAlt2<-intvarAlt2/sim.n
  intvarAlt3<-intvarAlt3/sim.n
  finvarAlt<-ifelse(identical(num.stage,2),intvarAlt2,intvarAlt3)
  pinfoAlt <- finvarAlt/intvarAlt
  if(identical(num.stage,3))pinfoAlt2 <- finvarAlt/intvarAlt2
  pstopAlt  <- pstopAlt/sim.n
  pstopEAlt  <- pstopEAlt/sim.n
  if(identical(num.stage,2)){
      n2int<-NA
      t2int<-NA
      aveE2<-NA
      pinfoNull2<-NA
      pinfoAlt2<-NA
  }


  names(t1int)<-""
  names(n1int)<-""
  names(t2int)<-""
  names(n2int)<-""
  names(eda)<-""
  names(etsl)<-""
  names(edaAlt)<-""
  names(etslAlt)<-""
  names(pinfoNull)<-""
  names(pinfoNull2)<-""
  names(pinfoAlt)<-""
  names(pinfoAlt2)<-""
  names(pstopNull)<-""
  names(pstopAlt)<-""
  names(pstopENull)<-""
  names(pstopEAlt)<-""
  names(es)<-""
  names(alphaExact)<-""
  names(alphaNorm)<-""
  names(powerExact)<-""
  names(powerNorm)<-""


  return(c(alphaExact=alphaExact,alphaNorm=alphaNorm,
           powerExact=powerExact,powerNorm=powerNorm,
           eda=eda,etsl=etsl,es=es,
           edaAlt=edaAlt,etslAlt=etslAlt,esAlt=esAlt,
           n1=n1int,n2=n2int,t1=t1int,t2=t2int,
           aveE=aveE,aveE2=aveE2,
           pinfoNull=pinfoNull, pinfoNull2=pinfoNull2,
           pinfoAlt=pinfoAlt,pinfoAlt2=pinfoAlt2,
           difIntFutL=difIntFutL,difIntSupH=difIntSupH,
           difIntFutH=difIntFutH,difIntSupL=difIntSupL,
           difFinSupL=difFinSupL,difFinFutH=difFinFutH,
           pstopNull=pstopNull,pstopAlt=pstopAlt,
           pstopENull=pstopENull,pstopEAlt=pstopEAlt))

}


ft<-function(ti,y,x,EW,pause,mda,pt1=NULL){
  ### only specify pt1 when computing second interim time
  if((!is.null(pt1)) && (pt1<mda))y[y>pt1]<-y[y>pt1]+pause

  if(ti<mda)pamt<-pause else pamt<-min(pause,max(mda+x-ti,0))

  y[y>ti]<-y[y>ti]+pamt
  tymin<-pmax(0,ti+pamt-y)
  return( sum(pmin(x,tymin))-EW )
}

