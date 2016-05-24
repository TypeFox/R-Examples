### Clustering

simulWeib.clst<-function(N,duration,lambda,rho,beta,rateC,df,min.futime)
{ 
  
  if(sum(df$cat_prop)!=1){
    cat("Error: proportions of patients do not sum to 1")
  }
  else if(N<nrow(df)){
    cat("Error: not enough number of patients")
  }
  else{
    expose<-NULL
    time<-NULL
    status<-NULL
    
    pats_grp<-rmultinom(n=1,size=N,prob=df$cat_prop)
    clst_id<-unlist(mapply(rep,df$cat_id,pats_grp))
    start<-rep(0,N)
    
    for(i in 1:length(pats_grp)){
      expose_i<-rbinom(n=sum(clst_id==df$cat_id[i]),size=1,prob=df$cat_exp.prop[i])
      v_i<-runif(n=sum(clst_id==df$cat_id[i]))
      Tlat<-(-log(v_i)/(lambda*exp(expose_i*beta)))^(1/rho)
      
      C<-rexp(n=sum(clst_id==df$cat_id[i]),rate=rateC)
      C<-pmin(C,rep(duration,length(C)))
      time_i<-pmin(Tlat,C)    
      status_i<-as.numeric(Tlat<=C)
      
      expose<-c(expose,expose_i)
      time<-c(time,time_i)
      status<-c(status,status_i)
    }
    if(min.futime==0){
      return(data.frame(id=1:length(time),start=start,stop=time,status=status,x=expose,
                        clst_id=clst_id))
    }
    else{
      return(data.frame(id=1:length(time),start=start,stop=time,status=status,x=expose,
                        clst_id=clst_id)[which(time>min.futime),])
    }  
  }
}

# modified version to generate time-dependent dataset with clustering
#' @export
tdSim.clst<-function(N,duration=24,lambda,rho=1,beta,rateC,df,
                     prop.fullexp=0,maxrelexptime=1,min.futime=0,min.postexp.futime=0){
  data<-simulWeib.clst(N,duration,lambda,rho,beta,rateC,df,min.futime)
  if(sum(df$cat_prop)==1 & N>=nrow(df)){
    if(prop.fullexp==0){
      data_tdexposed<-data[data$x==1,]
    }
    else{
      id_tdexposed<-sample(x = data[data$x==1,]$id,size = round(nrow(data[data$x==1,])*(1-prop.fullexp)))
      data_tdexposed<-data[data$id %in% id_tdexposed,]
    }
    data_tdexposed$t_exposed<-runif(nrow(data_tdexposed),0,data_tdexposed$stop*maxrelexptime)
    if(min.postexp.futime>0){
      if(sum(data_tdexposed$stop-data_tdexposed$t_exposed>min.postexp.futime) == 0){
        cat("Warning: no exposure left")
      }
      data_tdexposed<-data_tdexposed[data_tdexposed$stop-data_tdexposed$t_exposed>min.postexp.futime,]
    }
    new_data1<-data_tdexposed
    new_data2<-data_tdexposed
    new_data1$id<-data_tdexposed$id
    new_data1$start<-data_tdexposed$start
    new_data1$stop<-data_tdexposed$t_exposed
    new_data1$status<-0 
    new_data1$x<-0 
    new_data1$clst_id<-data_tdexposed$clst_id
    new_data2$id<-data_tdexposed$id
    new_data2$start<-data_tdexposed$t_exposed
    new_data2$stop<-data_tdexposed$stop
    new_data2$status<-data_tdexposed$status
    new_data2$x<-1
    new_data2$clst_id<-data_tdexposed$clst_id
    merged_tdexposed<-subset(na.omit(merge(new_data1,new_data2,all.x=TRUE,all.y=TRUE)))
	merged_tdexposed$t_exposed<-NULL
    full_data<-merge(merged_tdexposed,data[data$x==0,],all.x=TRUE,all.y=TRUE)
    return(full_data)
  }
}

# get.power function for clustering scenario
#' @export
getpower.clst<-function(nSim,N,duration=24,med.TTE.Control=24,rho=1,beta,med.TimeToCensor=14,df,type,scenario,
                        prop.fullexp=0,maxrelexptime=1,min.futime=0,min.postexp.futime=0,output.fn,simu.plot=FALSE) 
{
  lambda<-log(2)/med.TTE.Control
  rateC=log(2)/med.TimeToCensor
  #numsim=500
  res=matrix(0,nSim,8)
  colnames(res)=c("betahat","HR","signif","events",
                  "events_c","events_exp","medsurvt_c","medsurvt_exp")
  N.eff<-matrix(0,nSim,length(df$cat_id))
  N.effexp.p<-matrix(0,nSim,length(df$cat_id))
  alpha=.05
  if(simu.plot){
    set.seed(999)
    if(type == "fixed"){
      dat <- simulWeib.clst(N=N,duration=duration,lambda=lambda,rho=rho,beta=beta,rateC=rateC,
                            df=df,min.futime=min.futime)
    }
    else{
      dat <- tdSim.clst(N=N,duration=duration,lambda=lambda,rho=rho,beta=beta,rateC=rateC,
                        df=df,prop.fullexp=prop.fullexp,maxrelexptime=maxrelexptime,
                        min.futime=min.futime,min.postexp.futime=min.postexp.futime)
    }
    plot_simuData(dat)
  }
  set.seed(999)
  for(k in 1:nSim)
  {
    if(type == "fixed"){
      dat<-simulWeib.clst(N=N,duration=duration,lambda=lambda,rho=rho,beta=beta,rateC=rateC,
                          df=df,min.futime=min.futime)
    }
    else{
      dat<-tdSim.clst(N=N,duration=duration,lambda=lambda,rho=rho,beta=beta,rateC=rateC,
                      df=df,prop.fullexp=prop.fullexp,maxrelexptime=maxrelexptime,
                      min.futime=min.futime,min.postexp.futime=min.postexp.futime)  
    }
    fit <- coxph(Surv(start,stop, status) ~ factor(x)+cluster(clst_id), data=dat)
    sfit <- survfit(Surv(start,stop, status) ~ factor(x)+cluster(clst_id), data=dat)
    res[k,"betahat"] <- summary(fit)$coef[,"coef"]
    res[k,"HR"] <- summary(fit)$coef[,"exp(coef)"]
    res[k,"signif"] <- ifelse(summary(fit)$coef[,"Pr(>|z|)"]<alpha,1,0)
    res[k,"events"] <- sum(dat$status)
    res[k,"events_c"] <- summary(sfit)$table[1,'events']
    res[k,"events_exp"] <- summary(sfit)$table[2,'events']
    res[k,"medsurvt_c"] <- summary(sfit)$table[1,'median']
    res[k,"medsurvt_exp"] <- summary(sfit)$table[2,'median']
    for(j in 1:length(unique(dat$clst_id))){
      N.eff[k,j]<-length(unique(dat[dat$clst_id==df$cat_id[j],]$id))
      N.effexp.p[k,j]<-sum(dat[dat$clst_id==df$cat_id[j],]$x)/length(unique(dat[dat$clst_id==df$cat_id[j],]$id))
    }
  }
  df=data.frame(i_scenario=scenario,
                i_type=type,
                i_N=N,
                i_min.futime=min.futime,
                i_min.postexp.futime=min.postexp.futime,
                i_cat=df$cat_id,
                i_cat_prop=df$cat_prop,
                i_cat_exp.prop=df$cat_exp.prop,
                i_exp.prop=sum(df$cat_prop*df$cat_exp.prop),
                i_lambda=lambda,
                i_rho=rho,
                i_rateC=rateC,                       
                i_beta=beta,
                N_eff=colMeans(N.eff),
                N_effexp_p=colMeans(N.effexp.p),
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
