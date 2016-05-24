"print.OptimDes" <-
function(x,dig=3,all=FALSE,condPow=F,CMadj=F,...)
{
  alpha <- x$test[1]
  beta <- x$test[2]  
  shape0 <- x$test[3]
  scale0 <- x$test[4]
  shape1 <- x$test[5]
  scale1 <- x$test[6]
  xx <- x$test[7]
  r <- x$design[2]
  num.arm <- as.numeric(x$design[1])
  num.stage <- as.numeric(x$design[3])
  pause<-as.numeric(x$design['pause'])
  t1<-  as.numeric(x$stageTime["t1"])
  mtsl<- as.numeric(x$stageTime["MTSL"])
  mda<-as.numeric(x$stageTime['MDA'])
  etsl<-as.numeric(x$result['ETSL'])
  eda<-as.numeric(x$result['EDA'])
  es<-as.numeric(x$result['ES'])
  if(identical(num.stage,2)){
    n1<-as.numeric(x$n[1])
    ntot<-as.numeric(x$n[2])
    nvec<-c(n1,ntot)
    C1<-  x$boundary["C1L"]
  }else if(identical(num.stage,3)){
    n1<-as.numeric(x$n[1])
    n2<-as.numeric(x$n[2])
    ntot<-as.numeric(x$n[3])
    nvec<-c(n1,n2,ntot)
    t2 <- x$stageTime["t2"]
    C1L<-x$boundary["C1L"]
    C1U<-x$boundary["C1U"]
    C2L<-x$boundary["C2L"]
    C2U<-x$boundary["C1U"]
    C3<-x$boundary["C3"]
  }
  pstopNull<-x$result['pstopNull']
  pstopAlt<-x$result['pstopAlt']
  se<- x$se
  Stime<-x$single.stageTime
  exposure<-x$exposure[1]
  if(identical(num.stage,3)) exposure<-c(exposure,x$exposure[2])

  sf<- x$sf
  if(identical(sf,'futility'))sf<-'Futility only'
  else if(identical(sf,'OF'))sf<-'Futility plus Obrien-Fleming boundary'
  else if(identical(sf,'Pocock'))sf<-'Futility plus Pocock boundary'



  l0<-lambda(shape0,scale0,xx)
  l1<-lambda(shape1,scale1,xx)

  if(isTRUE(CMadj) && !identical(num.arm,1))    
    stop("The CM adjustment is for single-arm trials only.")
  if(isTRUE(CMadj) && identical(num.stage,3))  
    stop("The CM adjustment is for two-stage trials only.")  


  ### apply adjustment to exact binomial fixed design size
  ### per CM 2003
  if (CMadj){
        CMadjfac<-Stime["n0E"]/Stime["n0"]
        ntot<-CMadjfac*ntot
        n1<-CMadjfac*n1
        nvec<-c(n1,ntot)
        mda<-ifelse(t1<mda,CMadjfac*(mda-pause)+pause,CMadjfac*mda)
        mtsl<-mda+xx
        t1<-CMadjfac*t1
        if(t1>=mda){
           pamt<-min(pause,mda+x-t1)
        }else pamt<-pause

        etsl<-pstopNull*(t1+pamt)+(1-pstopNull)*mtsl
        eda<-ifelse(t1<mda,pstopNull*t1+(1-pstopNull)*mda,mda)
        es<-n1+(1-pstopNull)*(ntot-n1)
        exposure<-x$exposure[2]
  }else CMadjfac<-1
  
  cat("             Optimal Design Results \n\n")
  if(identical(num.arm,1) & identical(num.stage,2))
  cat("    One-Arm Two-Stage Study \n\n")
  if(identical(num.arm,1) & identical(num.stage,3))
  cat("    One-Arm Three-Stage Study \n\n")  
  if(identical(num.arm,2) & identical(num.stage,2))
  cat("    Two-Arm Two-Stage Study:", r/(1-r),": 1 randomization \n\n")
  if(identical(num.arm,2) & identical(num.stage,3))
  cat("    Two-Arm Three-Stage Study:", r/(1-r),": 1 randomization \n\n")

  cat("    Interim stopping rule:", sf,"\n\n")
  
  cat("    Pause in accrual before interim analyses:", pause,"\n\n")

  cat("    H0: S0=S1=",round(s(shape0,scale0,xx),dig),"H1: S1=",round(s(shape1,scale1,xx),dig),"\n\n")
  cat("    Type I error(1-sided upper):",alpha,"type II error:",beta,"\n")
  cat("    Event-free endpoint time:",xx," \n\n")
  
  cat("                  target:",x$target,"\n")
  prt1 <- c(eda,etsl,es)
  prt1<-round(prt1,dig)
  names(prt1) <- c("EDA","ETSL","           ES")
  print(prt1,...)

  cat("\n")
  cat("          Sample Size at Each Stage \n")
  prt2 <- ceiling(nvec)
  if(isTRUE(max(prt2)>floor(sum(unlist(x$accrual["m.init"]))))) 
          warning("Projected patient sample size is below the minimum requirement for normal approximation adjustment")
  
  if(identical(num.stage,2)){
  names(prt2) <- c("n1","          nmax")
  }else if(identical(num.stage,3)){
  names(prt2) <- c("n1","          n2","          nmax")
  }
  print(prt2,...)

  cat("\n")
  cat("          Study time at Each Stage \n")
  if(identical(num.stage,2)){
  prt3 <- round(c(t1,mda,mtsl),dig)
  names(prt3) <- c("t1","             MDA","             MTSL")
  }else if(identical(num.stage,3)){
  prt3 <- round(c(t1,t2,mda,mtsl),dig)
  names(prt3) <- c("t1","             t2",
                   "             MDA","             MTSL")
  }
  print(prt3,...)

  cat("\n")
  cat("    Projected patient exposure at interim analysis: ",
      round(exposure,2))
  cat("\n\n")
  cat("    Proportion of the total information at the interim analysis:\n")
  if(identical(num.arm,1) & identical(num.stage,2)){
  prt3b<-round(c((se[2]/se[1])^2,(se[4]/se[3])^2),3)
  names(prt3b)<-c("Under Null","    Under Alternative")  
  }
  if(identical(num.arm,1) & identical(num.stage,3)){
  prt3b<-round(c((se[3]/se[1])^2,(se[3]/se[2])^2,(se[6]/se[4])^2,(se[6]/se[5])^2),3)
  names(prt3b)<-c("Under Null Stage 1","    Under Null Stage 2","    Under Alternative Stage 1","    Under Alternative Stage 2")   
  }
  if(identical(num.arm,2) & identical(num.stage,2))
  {
    sig10 <- se[1]
    sig20 <- se[2]
    sig11 <- se[3]
    sig21 <- se[4]
    v10 <- (sig10^2)/((1-r)*(l0)^2)
    v11 <- (sig11^2)/(r*(l1)^2)
    v20 <- (sig20^2)/((1-r)*(l0)^2)
    v21 <- (sig21^2)/(r*(l1)^2)
    
    prt3b<-round(c((sig20/sig10)^2,(v20+v21)/(v10+v11)),3)
    names(prt3b)<-c("Under Null","    Under Alternative") 
  } 
  if(identical(num.arm,2) & identical(num.stage,3))
  {
    sig10 <- se[1]
    sig20 <- se[2]
    sig30 <- se[3]
    sig11 <- se[4]
    sig21 <- se[5]
    sig31 <- se[6]
    v10 <- (sig10^2)/((1-r)*(l0)^2)
    v11 <- (sig11^2)/(r*(l1)^2)
    v20 <- (sig20^2)/((1-r)*(l0)^2)
    v21 <- (sig21^2)/(r*(l1)^2)
    v30 <- (sig30^2)/((1-r)*(l0)^2)
    v31 <- (sig31^2)/(r*(l1)^2)
    
    prt3b<-round(c((sig30/sig10)^2,(sig30/sig20)^2,(v30+v31)/(v10+v11),(v30+v31)/(v20+v21)),3)
    names(prt3b)<-c("Under Null Stage 1","    Under Null Stage 2","    Under Alternative Stage 1","    Under Alternative Stage 2")
  }    
  

  print(prt3b,...)

  cat("\n\n")
  cat("          Hypothesis Test Boundaries \n")
  prt4 <- round(x$boundary,dig)
  if(identical(num.stage,2)){
  names(prt4) <- c("C1L","               C1U","               C2U")
  } else if(identical(num.stage,3)){
  names(prt4) <- c("C1L","               C1U","               C2L","               C2U","               C3U")
  }
  print(prt4,...)

  cat("\n\n")
  cat("    Approximate Rates Corresponding to Test Boundaries* \n\n")
  ### compute boundaries using the null and alternative SE's
  if(identical(num.arm,1) & identical(num.stage,2)) {
      se10<-x$se[1]/sqrt(ntot)  
      se20<-x$se[2]/sqrt(ntot)
      se11<-x$se[3]/sqrt(ntot)
      se21<-x$se[4]/sqrt(ntot)
      C1<- x$boundary[1]
      C2<- x$boundary[3]
      p10<-exp(log(l0)-C1*se10/l0)
      p10<-exp(-p10)
      p11<-exp(log(l0)-C1*se11/l1)
      p11<-exp(-p11)
      p1<-(p10+p11)/2

      p20<-exp(log(l0)-C2*se20/l0)
      p20<-exp(-p20)
      p21<-exp(log(l0)-C2*se21/l1)
      p21<-exp(-p21)
      p2<-(p20+p21)/2
      cat("          Event-free rate for C1L: ", 
          round(p1,dig),"\n")

      cat("          Event-free rate for C2U: ", 
          round(p2,dig),"\n\n")
      bndP<-c(p1,p2)  
      names(bndP)<-c('p1','p2')
  }
  
  if(identical(num.arm,1) & identical(num.stage,3)) {
      se10<-x$se[1]/sqrt(ntot)  
      se20<-x$se[2]/sqrt(ntot)
      se30<-x$se[3]/sqrt(ntot)
      se11<-x$se[4]/sqrt(ntot)
      se21<-x$se[5]/sqrt(ntot)
      se31<-x$se[6]/sqrt(ntot)
      C1<- x$boundary[1]
      C2<- x$boundary[3]
      C3<- x$boundary[5]
      p10<-exp(log(l0)-C1*se10/l0)
      p10<-exp(-p10)
      p11<-exp(log(l0)-C1*se11/l1)
      p11<-exp(-p11)
      p1<-(p10+p11)/2

      p20<-exp(log(l0)-C2*se20/l0)
      p20<-exp(-p20)
      p21<-exp(log(l0)-C2*se21/l1)
      p21<-exp(-p21)
      p2<-(p20+p21)/2
      
      p30<-exp(log(l0)-C3*se30/l0)
      p30<-exp(-p30)
      p31<-exp(log(l0)-C3*se31/l1)
      p31<-exp(-p31)
      p3<-(p30+p31)/2
      cat("          Event-free rate for C1L: ", 
          round(p1,dig),"\n")

      cat("          Event-free rate for C2L: ", 
          round(p2,dig),"\n\n")
          
      cat("          Event-free rate for C3U: ", 
          round(p3,dig),"\n\n")    
       bndP<-c(p1,p2,p3)  
       names(bndP)<-c('p1','p2','p3')
  }  
  
  
  
  if(identical(num.arm,2) & identical(num.stage,2)) {
      vnull1 <- (sig10^2)/((1-r)*r*(l0^2))
      vnull2 <- (sig20^2)/((1-r)*r*(l0^2))
      valt1  <- v10+v11
      valt2  <- v20+v21
      C1<- x$boundary[1] 
      C1U<- x$boundary[2]
      C2<- x$boundary[3]
      
      p10<-exp(log(l0)-C1*sqrt(vnull1/(ntot)))
      d10<-exp(-p10) - exp(-l0)
      p11<-exp(log(l0)-C1*sqrt(valt1/(ntot)))
      d11<-exp(-p11) - exp(-l0)
      d1<-(d10+d11)/2
      
      pu0<-exp(log(l0)-C1U*sqrt(vnull1/(ntot)))
      du0<-exp(-pu0) - exp(-l0)
      pu1<-exp(log(l0)-C1U*sqrt(valt1/(ntot)))
      du1<-exp(-pu1) - exp(-l0)
      d1u<-(du0+du1)/2
      
      p20<-exp(log(l0)-C2*sqrt(vnull2/(ntot)))
      d20<-exp(-p20) - exp(-l0)
      p21<-exp(log(l0)-C2*sqrt(valt2/(ntot)))
      d21<-exp(-p21) - exp(-l0)
      d2<-(d20+d21)/2
      cat("          Difference in Event-free rate for C1L: ", 
          round(d1,dig),"\n")
      if(!identical(sf,'Futility only')){
          cat("          Difference in Event-free rate for C1U: ", 
              round(d1u,dig),"\n")
      }
      cat("          Difference in Event-free rate for C2: ", 
          round(d2,dig),"\n\n")
      bndP<-c(d1=d1,d1u=d1u,d2=d2)
      names(bndP)<-c('d1','d1u','d2')
  }

  if(identical(num.arm,2) & identical(num.stage,3)) {
      vnull1 <- (sig10^2)/((1-r)*r*(l0^2))
      vnull2 <- (sig20^2)/((1-r)*r*(l0^2))
      vnull3 <- (sig30^2)/((1-r)*r*(l0^2))
      valt1  <- v10+v11
      valt2  <- v20+v21
      valt3  <- v30+v31
      C1<- x$boundary[1] 
      C1U<- x$boundary[2]
      C2<- x$boundary[3]
      C2U<- x$boundary[4]
      C3<- x$boundary[5] 
      
      p10<-exp(log(l0)-C1*sqrt(vnull1/(ntot)))
      d10<-exp(-p10) - exp(-l0)
      p11<-exp(log(l0)-C1*sqrt(valt1/(ntot)))
      d11<-exp(-p11) - exp(-l0)
      d1<-(d10+d11)/2
      
      pu0<-exp(log(l0)-C1U*sqrt(vnull1/(ntot)))
      du0<-exp(-pu0) - exp(-l0)
      pu1<-exp(log(l0)-C1U*sqrt(valt1/(ntot)))
      du1<-exp(-pu1) - exp(-l0)
      d1u<-(du0+du1)/2
      
      p20<-exp(log(l0)-C2*sqrt(vnull2/(ntot)))
      d20<-exp(-p20) - exp(-l0)
      p21<-exp(log(l0)-C2*sqrt(valt2/(ntot)))
      d21<-exp(-p21) - exp(-l0)
      d2<-(d20+d21)/2
      
      pu0<-exp(log(l0)-C2U*sqrt(vnull2/(ntot)))
      du0<-exp(-pu0) - exp(-l0)
      pu1<-exp(log(l0)-C2U*sqrt(valt2/(ntot)))
      du1<-exp(-pu1) - exp(-l0)
      d2u<-(du0+du1)/2 
      
      p30<-exp(log(l0)-C3*sqrt(vnull3/(ntot)))
      d30<-exp(-p30) - exp(-l0)
      p31<-exp(log(l0)-C3*sqrt(valt3/(ntot)))
      d31<-exp(-p31) - exp(-l0)
      d3<-(d30+d31)/2
      
      cat("          Difference in Event-free rate for C1L: ", 
          round(d1,dig),"\n")
      if(!identical(sf,'Futility only')){
          cat("          Difference in Event-free rate for C1U: ", 
              round(d1u,dig),"\n")
      }
      cat("          Difference in Event-free rate for C2L: ", 
          round(d2,dig),"\n\n")
      if(!identical(sf,'Futility only')){
          cat("          Difference in Event-free rate for C2U: ", 
              round(d2u,dig),"\n")
      }
      cat("          Difference in Event-free rate for C3U: ", 
          round(d3,dig),"\n\n")
      bndP<-c(d1=d1,d1u=d1u,d2=d2,d2u=d2u,d3=d3)
      names(bndP)<-c('d1','d1u','d2l','d2u','d3')

  }

  cat("\n\n")
  cat("    Probability of stopping at an interim analysis \n")
  cat("         Under the Null:",round(pstopNull,3),"\n")
  cat("         Under the Alternative:",round(pstopAlt,3),"\n")

  condpow1<-NULL
  condpow1u<-NULL
  if(condPow){    ##calculate conditional power for the final analysis based on the first interim for practical considerations
     cat("\n")

     # P(Z2>C2|Z1=C1,H1) based on bivariate normal distributions
     if(identical(num.arm,1) && identical(num.stage,2)){
          rho<-x$se[4]/x$se[3]     
          condpow1<-1-pnorm(C2,mean=x$u[2] + rho*(C1-x$u[1]),sd=sqrt(1-rho^2))
          cat("          Conditional power at interim test boundary under H1: ",
          round(condpow1,dig), "\n\n")
          }
          
          
     if(identical(num.arm,1) && identical(num.stage,3)){
          rho<-x$se[6]/x$se[4]
          condpow1<-1-pnorm(C3,mean=x$u[3] + rho*(C1-x$u[1]),sd=sqrt(1-rho^2))
          cat("          Conditional power at first interim test boundary under H1: ",
          round(condpow1,dig), "\n\n")
          }
     
     
     if(identical(num.arm,2) && identical(num.stage,2)){    
         rho<-sqrt((v20+v21)/(v10+v11))
         condpow1<-1-pnorm(C2,mean=x$u[2] + rho*(C1-x$u[1]),sd=sqrt(1-rho^2))
         condpow1u<-1-pnorm(C2,mean=x$u[2] + rho*(C1U-x$u[1]),sd=sqrt(1-rho^2))
          
         cat("          Conditional power at interim test boundary C1L under H1: ",
         round(condpow1,dig), "\n\n")
         cat("          Conditional power at interim test boundary C1U under H1: ",
         round(condpow1u,dig), "\n\n")
         }
     
     if(identical(num.arm,2) && identical(num.stage,3)){    
         rho<-sqrt((v30+v31)/(v10+v11))
         condpow1<-1-pnorm(C3,mean=x$u[3] + rho*(C1-x$u[1]),sd=sqrt(1-rho^2))
         condpow1u<-1-pnorm(C3,mean=x$u[3] + rho*(C1U-x$u[1]),sd=sqrt(1-rho^2))
          
         cat("          Conditional power at first interim test boundary C1L under H1: ",
         round(condpow1,dig), "\n\n")
         cat("          Conditional power at first interim test boundary C1U under H1: ",
         round(condpow1u,dig), "\n\n")
         }
  }

  cat("\n")
  comment <- ifelse(identical(num.arm,1),"(Exact binomial calculation)", "(Fisher exact calculation)")
  cat(paste("    Single-stage Design ",comment,sep=""), "\n")
  prt5 <- c(round(Stime[4]),round(Stime[c(5:6)],dig))
  names(prt5) <- c("      Single stage N","DA","          SL")
  print(prt5,...)
  cat("\n")
  comment <- "(Asymptotic normal calculation)"
  cat(paste("    Single-stage Design ",comment,sep=""), "\n")
  prt5 <- c(round(Stime[1]),round(Stime[c(2:3)],dig))
  names(prt5) <- c("      Single stage N","DA","          SL")
  print(prt5,...)

  if(isTRUE(all))
  {
    if(identical(num.stage,2)){ 
        all.info<-x$all.info
        all.info[,1]<-ceiling(all.info[,1])
        all.info<-all.info[,c(1:2,6:10)]
    }else if(identical(num.stage,3)){
        all.info<-x$all.info
        all.info<-all.info[,c(1:3,9:13)]
    }
    if(isTRUE(max(all.info[,1])>floor(sum(unlist(x$accrual["m.init"]))))) 
          warning("Projected patient sample size is below the minimum requirement for normal approximation adjustment")
    cat("\n")
    cat("   Design Parameters for Each n \n")
    print(round(all.info,dig))
    if(CMadj)cat('\n Note: CM adjustment not included in the design parameters')
  }

  cat("\n\n *Note: Rates corresponding to test boundaries are a function \n",
      "of the non-parametric SE computed at the time of the analyses. \n", 
      "The approximate rates are based on the asymptotic SE computed  \n",
      "under the null and alternative hypotheses.\n")

  if(pause>0){
  cat(paste("\n *Note: Interim analysis time(s) are at the beginning of the\n",
            " accrual pause.  Information/exposure are computed at the end\n",
            " of the pause.\n",sep=''))
  }

  if(CMadj){
     cat("\n Note:  All sample sizes and times are adjusted by the exact ",
            "binomial correction \n","factor: ",Stime["n0E"],"/",
            Stime["n0"],"\n",sep="")
  }
  
  ret<-list(pinfo=prt3b,condpow=c(condpow1,condpow1u),
            bndP=bndP,CMadjfac=CMadjfac)
  return(invisible(ret))
}  
