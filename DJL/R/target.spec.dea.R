target.spec.dea <-
function(xdata,ydata,date=NULL,t=NULL,dt=NULL,dmu,et="c",alpha=NULL,beta=NULL,wv,rts,sg="ssm",ftype="d",ncv=NULL,env=NULL){
  
  # Initial checks
  if(is.null(date) && sg!="ssm"){stop('sg must be "ssm" when date is null.')}
  if(is.na(match(rts,c("crs","vrs","irs","drs")))){stop('rts must be "crs", "vrs", "irs", or "drs".')}
  if(is.na(match(sg,c("ssm","max","min")))){stop('sg must be "ssm", "max", or "min".')}
  if(is.na(match(ftype,c("d","s")))){stop('ftype must be either "d" or "s".')}
  if(dmu>nrow(xdata)||dmu<1){stop('dmu must indicate one of data points in the set.')}
  if(!xor(is.null(alpha),is.null(beta))){stop('Either alpha or beta must be defined.')}
  if(is.null(date)||is.null(t)||is.null(dt)) mtype<-"sidea" else mtype<-"tidea"
  if(mtype=="tidea" && t<=min(date)){stop('t is earlier than dataset.')}
  if(mtype=="tidea" && max(date)<t){stop('t is later than dataset.')}

  # Estimation of the orientation
  xdata<-as.matrix(xdata);ydata<-as.matrix(ydata);if(mtype=="tidea")date<-as.matrix(date) # format input data as matrix
  m<-ncol(xdata); s<-ncol(ydata)
  if(is.null(alpha) && !is.null(beta)){alpha<-matrix(rep(NA,m),nrow=1,ncol=m);orientation<-"i"}
  if(!is.null(alpha) && is.null(beta)){beta<-matrix(rep(NA,s),nrow=1,ncol=s);orientation<-"o"}
  if(orientation=="i" && ncol(as.matrix(wv))!=m){stop('wv must have the same number of column with xdata.')}
  if(orientation=="o" && ncol(as.matrix(wv))!=s){stop('wv must have the same number of column with ydata.')}
  if(!is.null(env)) env<-as.matrix(env)
  
  # Construction of the PPS
  if(mtype=="tidea"){
    # Calc RoCs
    grapejelly<-roc.dea(xdata,ydata,date,t,rts,orientation,sg,ftype)
    if(et=="c"){et<-grapejelly$eff_t[dmu,]}
    
    # Future reference set
    soa<-which(grapejelly$roc_local>0)
    d_f<-as.matrix(date[soa,])
    if(orientation=="i"){x_f<-as.matrix(xdata[soa,])*((1/grapejelly$roc_local[soa,])^dt);y_f<-as.matrix(ydata[soa,])}
    if(orientation=="o"){x_f<-as.matrix(xdata[soa,]);y_f<-as.matrix(ydata[soa,])*grapejelly$roc_local[soa,]^dt}
    if(!is.null(env)){env_f<-as.matrix(env[soa,])}
    
    # Experiment set
    d_e<-rbind(d_f,date[dmu,])
    if(orientation=="i"){x_e<-rbind(x_f,xdata[dmu,]);y_e<-rbind(y_f,beta)}
    if(orientation=="o"){x_e<-rbind(x_f,alpha);y_e<-rbind(y_f,ydata[dmu,])}
    if(!is.null(env)){env_e<-rbind(env_f,env[dmu,])}
  }else{
    # Calc efficiency
    grapejelly<-dm.dea(xdata,ydata,rts,orientation)
    if(et=="c"){et<-grapejelly$eff[dmu,]}
    
    # Experiment set
    soa<-which(round(grapejelly$eff,8)==1)
    if(orientation=="i"){x_e<-rbind(as.matrix(xdata[soa,]),xdata[dmu,]);y_e<-rbind(as.matrix(ydata[soa,]),beta)}
    if(orientation=="o"){x_e<-rbind(as.matrix(xdata[soa,]),alpha);y_e<-rbind(as.matrix(ydata[soa,]),ydata[dmu,])}
    if(!is.null(env)){env_e<-rbind(as.matrix(env[soa,]),env[dmu,])}
  }
  
  # Parameters
  n<-nrow(x_e)
  if(is.null(ncv)) ncv<-matrix(c(0),ncol=m+s) else ncv<-as.matrix(ncv)
  
  # Feasibility check
  if(orientation=="i"){
    x_l<-rbind(as.matrix(x_e[1:(n-1),]),xdata[dmu,]*et)
    y_l<-rbind(as.matrix(y_e[1:(n-1),]),ydata[dmu,])
    bound<-ydata[dmu,]*(dm.dea(x_l,y_l,rts,"o")$eff[n])
    if(sum(beta>bound)>0){stop(paste0('Beta(',paste(beta,collapse=", "),') is greater than feasible bound(',paste(round(bound,4),collapse=", "),').'))}
  }
  if(orientation=="o"){
    x_l<-rbind(as.matrix(x_e[1:(n-1),]),xdata[dmu,])
    y_l<-rbind(as.matrix(y_e[1:(n-1),]),ydata[dmu,]*et)
    bound<-xdata[dmu,]*(dm.dea(x_l,y_l,rts,"i")$eff[n])
    if(sum(alpha<bound)>0){stop(paste0('Alpha(',paste(alpha,collapse=", "),') is smaller than feasible bound(',paste(round(bound,4),collapse=", "),').'))}
  }
  
  # Data frames
  lambda<-matrix(rep(NA,n),nrow=1,ncol=n)
  xslack<-matrix(rep(NA,m),nrow=1,ncol=m) 
  yslack<-matrix(rep(NA,s),nrow=1,ncol=s) 
  
  # Declare LP
  if(orientation=="i"){lp.idea<-make.lp(0,n+m+m+s)} # lambda+alpha+xslack+yslack
  if(orientation=="o"){lp.idea<-make.lp(0,n+s+m+s)} # lambda+beta+xslack+yslack
  
  # Set objective
  if(orientation=="i"){set.objfn(lp.idea,c(wv),indices=c((n+1):(n+m)))}
  if(orientation=="o"){set.objfn(lp.idea,c(-wv),indices=c((n+1):(n+s)))}
  
  # RTS
  if(rts=="vrs"){add.constraint(lp.idea,c(rep(1,n)),indices=c(1:n),"=",1)}
  if(rts=="crs"){set.constr.type(lp.idea,0,1)}
  if(rts=="irs"){add.constraint(lp.idea,c(rep(1,n)),indices=c(1:n),">=",1)}
  if(rts=="drs"){add.constraint(lp.idea,c(rep(1,n)),indices=c(1:n),"<=",1)}
  
  # Input constraints
  for(i in 1:m){
    if(orientation=="i"){add.constraint(lp.idea,c(x_e[,i],-et^(1-ncv[1,i]),1),indices=c(1:n,n+i,n+m+i),"=",0)}
    if(orientation=="o"){add.constraint(lp.idea,c(x_e[,i],1),indices=c(1:n,n+s+i),"=",alpha[1,i])}
  }
  
  # Output constraints
  for(r in 1:s){
    if(orientation=="i"){add.constraint(lp.idea,c(y_e[,r],-1),indices=c(1:n,n+m+m+r),"=",beta[1,r])}
    if(orientation=="o"){add.constraint(lp.idea,c(y_e[,r],-et^(1-ncv[1,i]),-1),indices=c(1:n,n+r,n+s+m+r),"=",0)}
  }
  
  # External NDF
  if(!is.null(env)){for(j in 1:n){if(env_e[j,1]<env_e[n,1]){add.constraint(lp.idea,c(1),indices=c(j),"=",0)}}}
  
  # Bounds
  if(orientation=="i"){set.bounds(lp.idea,lower=c(rep(0,n+m+m+s)))
    set.bounds(lp.idea,upper=c(rep(Inf,n),xdata[dmu,],rep(Inf,m+s)))}
  if(orientation=="o"){set.bounds(lp.idea,lower=c(rep(0,n),ydata[dmu,],rep(0,m+s)))}
  
  # Solve
  solve.lpExtPtr(lp.idea)
  
  # Get results
  temp.p<-get.variables(lp.idea)
  lambda[1,]<-temp.p[1:n]
  if(orientation=="i"){alpha[1,]<-temp.p[(n+1):(n+m)];xslack[1,]<-temp.p[(n+m+1):(n+m+m)];yslack[1,]<-temp.p[(n+m+m+1):(n+m+m+s)]}
  if(orientation=="o"){beta[1,]<-temp.p[(n+1):(n+s)];xslack[1,]<-temp.p[(n+s+1):(n+s+m)];yslack[1,]<-temp.p[(n+s+m+1):(n+s+m+s)]}
  
  # Stage II
  if(exists("sg")){
    # Link previous solutions
    if(orientation=="i"){for(i in 1:m){add.constraint(lp.idea,c(1),indices=c(n+i),"=",alpha[1,i])}}
    if(orientation=="o"){for(r in 1:s){add.constraint(lp.idea,c(1),indices=c(n+r),"=",beta[1,r])}}
    
    # date sum
    if(sg=="max"){set.objfn(lp.idea,c(-d_e[1:n,]),indices=c(1:n))}
    if(sg=="min"){set.objfn(lp.idea,c(d_e[1:n,]),indices=c(1:n))}
    
    # slack sum max
    if(sg=="ssm"){
      if(orientation=="i"){set.objfn(lp.idea,c(rep(-1,m+s)),indices=c((n+m+1):(n+m+m+s)))}
      if(orientation=="o"){set.objfn(lp.idea,c(rep(-1,m+s)),indices=c((n+s+1):(n+s+m+s)))}
    }
    
    # solve
    solve.lpExtPtr(lp.idea)
    
    # get results
    temp.d<-get.variables(lp.idea)
    lambda[1,]<-temp.d[1:n]
    if(orientation=="i"){alpha[1,]<-temp.d[(n+1):(n+m)];xslack[1,]<-temp.d[(n+m+1):(n+m+m)];yslack[1,]<-temp.d[(n+m+m+1):(n+m+m+s)]}
    if(orientation=="o"){beta[1,]<-temp.d[(n+1):(n+s)];xslack[1,]<-temp.d[(n+s+1):(n+s+m)];yslack[1,]<-temp.d[(n+s+m+1):(n+s+m+s)]}
  }
  results<-list(alpha=alpha,beta=beta,lambda=lambda,xslack=xslack,yslack=yslack)
  return(results)    
}
