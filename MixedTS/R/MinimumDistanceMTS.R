MinDist.VG<-function(dataset,param0,N=50,MixingDens="Gamma"){
#   # The qmle function works in a similar way of the mle function
#   
  env<-new.env()
  env$MixingDens<-MixingDens
  env$data<-dataset
  env$N<-N
  
  MinDis.VG<-function(par,dataset,N){
    mu0<-par[1]
    mu<-par[2]
    sig<-par[3]
    a<-par[4]
    alpha<-2
    lambda_p<-1
    lambda_m<-1
    # #     dens<-na.omit(dMixedTS(env$data,mu0,
    # #                            mu,sig,a,alpha,lambda_p,lambda_m))
    # #     #     dens[is.na(dens)]<-1 # in this way we remove 
    #     
    cond<-is.finite(dataset)
    sol<-hist(as.numeric(dataset[cond]),nclass=N,plot = FALSE)
    xstep<-sol$mids
    densEmp<-sol$density
    densTheo<-dMixedTS(xstep,mu0,
                       mu,sig,a,alpha,lambda_p,lambda_m)
    cond2<-is.finite(densTheo)
    
    sum((densEmp[cond2]-densTheo[cond2])^2)
    #     #     dens[is.na(dens)]<-1 # in this way we remove 
    #     -sum(log(dens)[is.finite(log(dens))])   
    
  }
  if(MixingDens=="Gamma"){
    #  ui<-rbind(c(1, -1, 0, 0),c(1, 1, 0, 0),c(1, 0, 0, 0),c(0, 0, 1, 0))
    #   mu0,mu,sig,a,alpha,lambda_p,lambda_m
    ui<-rbind(c(0, 0, 1, 0 ),
              c(0, 0, 0, 1))
    ci<-c(10^(-6), 10^(-6))
    #  ci<-c(0,0,0,10^(-6))
    # We have to insert the parameters restrictions considered in the paper    
  }
  
  lengpar<-length(param0)
  paramLev<-NA*c(1:length(lengpar))
#   
  env$lengpar<-lengpar
  firs.prob<-tryCatch(constrOptim(theta=param0,
                                  f=MinDis.VG,grad=NULL,ui=ui,ci=ci,dataset=dataset,N=N),
                      error=function(theta){NULL})
  return(firs.prob)
#   
 }


MinDist.MixedTS<-function(data,param0, method="L-BFGS-B", fixed.param=NULL,
         lower.param=NULL,
         upper.param=NULL,
         MixingDens="Gamma",N=50){
  # The qmle function works in a similar way of the mle function
  
  env<-new.env()
  env$MixingDens<-MixingDens
  env$data<-data
  env$N<-N
  if(!is.null(fixed.param)){
    if(fixed.param["alpha"]==2){
      resDisVG<-MinDist.VG(dataset=data,param0=c(param0["mu0"],param0["mu"],
                                       param0["sig"],param0["a"]),N=N,
                           MixingDens=MixingDens)
      return(resDisVG)
    }
  }
  
  env$cond<-is.finite(env$data)
  env$sol<-hist(as.numeric(env$data[env$cond]),nclass=env$N,plot = FALSE)
  
  if(MixingDens=="Gamma"){
    #  ui<-rbind(c(1, -1, 0, 0),c(1, 1, 0, 0),c(1, 0, 0, 0),c(0, 0, 1, 0))
    #   mu0,mu,sig,a,alpha,lambda_p,lambda_m
    ui<-rbind(c(0, 0, 1, 0, 0, 0, 0),
              c(0, 0, 0, 1, 0, 0, 0),
              c(0, 0, 0, 0, 1, 0, 0),
              c(0, 0, 0, 0,-1, 0, 0),
              c(0, 0, 0, 0, 0, 1, 0),
              c(0, 0, 0, 0, 0, 0, 1))
    ci<-c(10^(-6), 10^(-6), 10^(-6), -(2-10^(-6)), 10^(-6), 10^(-6))
    #  ci<-c(0,0,0,10^(-6))
    # We have to insert the parameters restrictions considered in the paper    
  }
  if(!is.null(lower.param)){
    lower.con<-matrix(0,length(lower.param),length(param0))
    rownames(lower.con)<-names(lower.param)
    colnames(lower.con)<-names(param0)
    numb.lower<-length(lower.param)
    lower.con[names(lower.param),names(lower.param)]<-1*diag(numb.lower)
    dummy.lower.names<-paste0(names(lower.param),".lower")
    rownames(lower.con)<-dummy.lower.names
    names(lower.param)<-dummy.lower.names
    ui<-rbind(ui,lower.con)
    ci<-c(ci,lower.param)
    #idx.lower.param<-match(names(lower.param),names(param0))
  }
  if(!is.null(upper.param)){
    upper.con<-matrix(0,length(upper.param),length(param0))
    rownames(upper.con)<-names(upper.param)
    colnames(upper.con)<-names(param0)
    numb.upper<-length(upper.param)
    upper.con[names(upper.param),names(upper.param)]<--1*diag(numb.upper)
    dummy.upper.names<-paste0(names(upper.param),".upper")
    rownames(upper.con)<-dummy.upper.names
    names(upper.param)<-dummy.upper.names
    ui<-rbind(ui,upper.con)
    ci<-c(ci,-upper.param)
  }
  if(!is.null(fixed.param)){
    names.fixed<-names(fixed.param)
    numb.fixed<-length(fixed.param)
    fixed.con<-matrix(0,length(fixed.param),length(param0))
    rownames(fixed.con)<-names(fixed.param)
    colnames(fixed.con)<-names(param0)
    fixed.con.bis<-fixed.con
    fixed.con[names(fixed.param),names(fixed.param)]<--1*diag(numb.fixed)
    fixed.con.bis[names(fixed.param),names(fixed.param)]<-1*diag(numb.fixed)
    dummy.fixed.names<-paste0(names(fixed.param),".fixed.u")
    dummy.fixed.bis.names<-paste0(names(fixed.param),".fixed.l")
    rownames(fixed.con)<-dummy.fixed.names
    rownames(fixed.con.bis)<-dummy.fixed.bis.names
    names(fixed.param)<-dummy.fixed.names
    ui<-rbind(ui,fixed.con,fixed.con.bis)
    ci<-c(ci,-fixed.param-10^-6,fixed.param-10^-6)
    #ci<-c(ci,-fixed.param,fixed.param)
  }  
  
  lengpar<-length(param0)
  paramLev<-NA*c(1:length(lengpar))
  
  env$lengpar<-lengpar
  time<-system.time(
  firs.prob<-tryCatch(constrOptim(theta=param0,
                                  f=MinDis,grad=NULL,ui=ui,ci=ci,env=env,method=method),
                      error=function(theta){NULL})
  )
  if(!is.null(firs.prob)){
    paramLev<-firs.prob$par
    names(paramLev)<-names(param0)
    if(!is.null(fixed.param)){
      paramLev[names.fixed]<-fixed.param
      names(paramLev)<-names(param0)
    }
  }else{warning("the start value for levy measure is outside of the admissible region")}
  
  
  results<-list(estLevpar=paramLev,
                MeasureDist=firs.prob$value,time=time)
  
  return(results)
}

MinDis<-function(par,env){
       mu0<-par[1]
       mu<-par[2]
       sig<-par[3]
       a<-par[4]
       alpha<-par[5]
       lambda_p<-par[6]
       lambda_m<-par[7]
  # #     dens<-na.omit(dMixedTS(env$data,mu0,
  # #                            mu,sig,a,alpha,lambda_p,lambda_m))
  # #     #     dens[is.na(dens)]<-1 # in this way we remove 
  #     
  #  xstep<-env$sol$mids
#  densEmp<-env$sol$density
    densTheo<-dMixedTS(env$sol$mids,mu0,
                             mu,sig,a,alpha,lambda_p,lambda_m,N=2^7)
#    cond2<-is.finite(densTheo)
  
  sum((env$sol$density[is.finite(densTheo)]-densTheo[is.finite(densTheo)])^2)
  #     #     dens[is.na(dens)]<-1 # in this way we remove 
  #     -sum(log(dens)[is.finite(log(dens))])   
  
}