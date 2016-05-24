#-----------------------------------------
# File name: ll.gd.R
#-----------------------------------------
ll.gd<-function(counts,r){
  
  # gene-dosage model
  C1=c(2,2,1,2,1,1,1,2,1,0,1,0,1,0,0)

  # Mating type variables
  MT1=c(1,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
  MT2=c(0,1,1,1,1,0,0,0,0,0,0,0,0,0,0)
  MT3=c(0,0,0,0,0,1,1,0,0,0,0,0,0,0,0)
  MT4=c(0,0,0,0,0,0,0,1,1,1,0,0,0,0,0)
  MT5=c(0,0,0,0,0,0,0,0,0,0,1,1,1,1,0)
  MT6=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,1)
  
  # offset for transmission probability -- TRD
  offset.trd=c(log(1),
               log(r),log(1-r),log(r),log(1-r),
               log(1),log(1),
               log(r^2),log(2*r*(1-r)),log((1-r)^2),
               log(r),log(1-r),log(r),log(1-r),
               log(1))
  
  # offset for M=F=C=1
  offset.111= c(rep(0,8),log(2),rep(0,6))
    
  # Fit log-linear model without TRD adjustment
  mod.0=glm(counts~MT1+MT2+MT3+MT4+MT5, 
            family=poisson,offset=offset.111)
  
  mod.wb=glm(counts~MT1+MT2+MT3+MT4+MT5+C1,
             family=poisson,offset=offset.111) 
  
  # Fit log-linear model with TRD adjustment
  mod.trd.0=glm(counts~MT1+MT2+MT3+MT4+MT5, 
                family=poisson,
                offset=offset.trd)
  
  mod.trd.wb=glm(counts~MT1+MT2+MT3+MT4+MT5+C1, 
                 family=poisson, 
                 offset=offset.trd)
    
  # Acquire deviance for the 2 models  
  dev.wb=anova(mod.0,mod.wb)$Deviance[2]
  dev.trd.wb=anova(mod.trd.0,mod.trd.wb)$Deviance[2]
  
  stat=list(
    summary(mod.wb)$coeff[7,],
    summary(mod.trd.wb)$coeff[7,],
    dev.wb,
    dev.trd.wb
  )
  
  names(stat)=c('M1','M2','DevM1','DevM2')
  stat
  
}