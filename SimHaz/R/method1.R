### method 1

simulWeib <- function(N, duration, lambda, rho, beta, rateC,exp.prop, min.futime)
{ 
  # covariate --> N Bernoulli trials
  expose<-rbinom(n=N,size=1,prob=exp.prop)
  
  # Weibull latent event times
  v <- runif(n=N)
  Tlat <- (- log(v) / (lambda * exp(expose * beta)))^(1 / rho)
  
  # censoring times
  C <- rexp(n=N, rate=rateC)
  C=pmin(C,rep(duration,length(C)))
  
  # follow-up times and event indicators
  time <- pmin(Tlat, C)
  
  status <- as.numeric(Tlat <= C)
  start = rep(0,length(time)) #all start at 0
  
  if(min.futime==0){
    return(data.frame(id=1:N,start=start,stop=time,status=status,x=expose))
  }
  else{
    return(data.frame(id=1:N,start=start,stop=time,status=status,x=expose)[which(time>min.futime),])
  }
}

# regular version to generate time-dependent dataset
# fullyexp.p: fully exposed proportion, the default value is 0, can take values in [0, 1)
# maxrelexp.t: maximum relative exposuret time, the default value is 1, can take values in (0, 1]
# min.postexp.fut: minimum post-exposure follow-up time
#' @export
tdSim.method1<-function(N, duration=24,lambda, rho=1, beta, rateC,exp.prop,
                        prop.fullexp=0,maxrelexptime=1,min.futime=0, min.postexp.futime=0){
  data <- simulWeib(N, duration,lambda, rho, beta, rateC,exp.prop,min.futime)
  if(prop.fullexp==0){
    data_tdexposed<-data[data$x==1,] #####add comment
  }
  else{
    id_tdexposed<-sample(x = data[data$x==1,]$id,size = round(nrow(data[data$x==1,])*(1-prop.fullexp)))
    data_tdexposed<-data[data$id %in% id_tdexposed,]
  }
  
  data_tdexposed$t_exposed<-runif(nrow(data_tdexposed),0,data_tdexposed$stop*maxrelexptime)
  
  if(min.postexp.futime>0){
    if(sum(data_tdexposed$stop-data_tdexposed$t_exposed>min.postexp.futime) == 0){
      print('no exposure left')
    }
    data_tdexposed <- data_tdexposed[data_tdexposed$stop-data_tdexposed$t_exposed>min.postexp.futime,]
  }
  new_data1<-data_tdexposed
  new_data2<-data_tdexposed
  new_data1$id<-data_tdexposed$id
  new_data1$start<-data_tdexposed$start
  new_data1$stop<-data_tdexposed$t_exposed
  new_data1$status<-0
  new_data1$x<-0
  new_data2$id<-data_tdexposed$id
  new_data2$start<-data_tdexposed$t_exposed
  new_data2$stop<-data_tdexposed$stop 
  new_data2$x<-1
  new_data2$status<-data_tdexposed$status      
  merged_tdexposed<-subset(na.omit(merge(new_data1,new_data2,all.x=TRUE,all.y=TRUE))) 
  merged_tdexposed$t_exposed<-NULL
  full_data<-merge(merged_tdexposed,data[data$x==0,],all.x=TRUE,all.y=TRUE)
  return(full_data)
}

#' @export
getpower.method1<-function(nSim, N,duration=24,med.TTE.Control=24,rho=1,med.TimeToCensor=14,beta,exp.prop,type,scenario,
                           prop.fullexp=0,maxrelexptime=1,min.futime=0,min.postexp.futime=0,output.fn,simu.plot=FALSE) 
{ 
  lambda=log(2)/med.TTE.Control
  rateC=log(2)/med.TimeToCensor
  #numsim=500
  res=matrix(0,nSim,10)
  colnames(res)=c("N.eff","N.effexp.p","betahat","HR","signif","events",
                  "events_c","events_exp","medsurvt_c","medsurvt_exp")
  alpha=.05
  if(simu.plot){
    set.seed(999)
    if(type == "fixed"){
      dat <- simulWeib(N=N, duration=duration,lambda=lambda, rho=rho, beta=beta, rateC=rateC,
                       exp.prop=exp.prop,min.futime=min.futime)
    }
    else{
      dat <- tdSim.method1(N=N, duration=duration,lambda=lambda, rho=rho, beta=beta, rateC=rateC,
                           exp.prop=exp.prop,prop.fullexp=prop.fullexp,maxrelexptime=maxrelexptime,
                           min.futime=min.futime,min.postexp.futime=min.postexp.futime)
    }
    plot_simuData(dat)
  }
  set.seed(999)
  for(k in 1:nSim)
  {
    if(type == "fixed"){
      dat <- simulWeib(N=N, duration=duration,lambda=lambda, rho=rho, beta=beta, rateC=rateC,
                       exp.prop=exp.prop,min.futime=min.futime)
    }
    else{
      dat <- tdSim.method1(N=N, duration=duration,lambda=lambda, rho=rho, beta=beta, rateC=rateC,
                           exp.prop=exp.prop,prop.fullexp=prop.fullexp,maxrelexptime=maxrelexptime,
                           min.futime=min.futime,min.postexp.futime=min.postexp.futime)     
    }
    fit <- coxph(Surv(start,stop, status) ~ factor(x), data=dat)
    sfit <- survfit(Surv(start,stop, status) ~ factor(x), data=dat)
    res[k,"N.eff"] <- length(unique(dat$id))
    res[k,"N.effexp.p"] <- sum(dat$x)/length(unique(dat$id))
    res[k,"betahat"] <- summary(fit)$coef[,"coef"]
    res[k,"HR"] <- summary(fit)$coef[,"exp(coef)"]
    res[k,"signif"] <- ifelse(summary(fit)$coef[,"Pr(>|z|)"]<alpha,1,0)
    res[k,"events"] <- sum(dat$status)
    res[k,"events_c"] <- summary(sfit)$table[1,'events']
    res[k,"events_exp"] <- summary(sfit)$table[2,'events']
    res[k,"medsurvt_c"] <- summary(sfit)$table[1,'median']
    res[k,"medsurvt_exp"] <- summary(sfit)$table[2,'median']
  }
  df=data.frame(i_scenario=scenario,
                i_type=type,
                i_N=N,
                i_min.futime=min.futime,
                i_min.postexp.futime=min.postexp.futime,
                i_exp.prop=exp.prop,
                i_lambda=lambda,
                i_rho=rho,
                i_rateC=rateC,                       
                i_beta=beta,
                N_eff=mean(res[,"N.eff"]),
                N_effexp_p=mean(res[,"N.effexp.p"]),
                bhat=mean(res[,"betahat"]),
                HR=mean(res[,"HR"]),                     
                d=mean(res[,"events"]),
                d_c=mean(res[,"events_c"]),
                d_exp=mean(res[,"events_exp"]),
                mst_c=mean(na.omit(res[,"medsurvt_c"])),
                mst_exp=mean(na.omit(res[,"medsurvt_exp"])),
                pow=mean(res[,"signif"])
  )
  if(file.exists(output.fn)){
    write.table(df,file=output.fn,row.names=FALSE,col.names=FALSE,append=TRUE,sep=",")
  }
  else{
    write.table(df,file=output.fn,row.names=FALSE,col.names=TRUE,sep=",")
  }
  return(df)
}
