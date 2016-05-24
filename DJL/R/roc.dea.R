roc.dea <-
function(xdata,ydata,date,t,rts,orientation,sg="ssm",ftype="d",ncv=NULL,env=NULL){

  # Initial checks
  if(is.na(match(rts,c("crs","vrs","irs","drs")))){stop('rts must be "crs", "vrs", "irs", or "drs".')}
  if(is.na(match(orientation,c("i","o")))){stop('orientation must be either "i" or "o".')}
  if(is.na(match(sg,c("ssm","max","min")))){stop('sg must be "ssm", "max", or "min".')}
  if(is.na(match(ftype,c("d","s")))){stop('ftype must be either "d" or "s".')}
  if(t<=min(date)){stop('t is earlier than dataset.')}
  if(max(date)<t){stop('t is later than dataset.')}
  
  # Parameters
  xdata<-as.matrix(xdata);ydata<-as.matrix(ydata);date<-as.matrix(date) # format input data as matrix
  n<-nrow(xdata); m<-ncol(xdata); s<-ncol(ydata)
  o<-matrix(c(1:n),ncol=1) # original data order
  
  # Sort data ascending order
  x<-matrix(c(xdata[order(date),]),ncol=m);colnames(x)<-colnames(xdata)
  y<-matrix(c(ydata[order(date),]),ncol=s);colnames(y)<-colnames(ydata)
  d<-matrix(c(date[order(date),]),ncol=1);colnames(d)<-colnames(date)
  o<-matrix(c(o[order(date),]),ncol=1)
  
  # Data frames
  eff_r<-array(NA,c(n,1))
  eff_t<-array(NA,c(n,1))
  lambda<-array(NA,c(n,n))
  ed<-array(NA,c(n,1))
  sl<-array(NA,c(n,1))
  roc<-array(NA,c(n,1))
  local_roc<-array(NA,c(n,1))
  
  # Subset index
  till<-function(x,y){
    t<-0
    while(x[t+1]<=y&&t<nrow(x)){t<-t+1}
    return(t)
  }
  r<-till(unique(d),t)

  # DEA internal function
  dm.dea.internal<-function(xdata,ydata,rts,orientation,se=0,sg,date,ncv,env,a,z){
    
    # Load library
    # library(lpSolveAPI)  
    
    # Parameters
    n<-nrow(xdata); m<-ncol(xdata); s<-ncol(ydata)
    if(is.null(ncv)){ncv<-matrix(c(0),ncol=m+s)}
    
    # Data frames
    results.efficiency<-matrix(rep(NA,n),nrow=n,ncol=1)
    results.lambda<-matrix(rep(NA,n^2),nrow=n,ncol=n)
    results.xslack<-matrix(rep(NA,n*m),nrow=n,ncol=m) 
    results.yslack<-matrix(rep(NA,n*s),nrow=n,ncol=s) 
    results.vweight<-matrix(rep(NA,n*m),nrow=n,ncol=m) 
    results.uweight<-matrix(rep(NA,n*s),nrow=n,ncol=s) 
    results.w<-matrix(rep(NA,n),nrow=n,ncol=1)
    
    for (k in a:z){
      # Declare LP
      lp.dea<-make.lp(0,n+1+m+s) # lambda+efficiency+xslack+yslack
      
      # Set objective
      if(orientation=="o"){set.objfn(lp.dea,c(-1),indices=c(n+1))}
      if(orientation=="i"){set.objfn(lp.dea,c(1),indices=c(n+1))}
      
      # RTS
      if(rts=="vrs"){add.constraint(lp.dea,c(rep(1,n)),indices=c(1:n),"=",1)}
      if(rts=="crs"){set.constr.type(lp.dea,0,1)}
      if(rts=="irs"){add.constraint(lp.dea,c(rep(1,n)),indices=c(1:n),">=",1)}
      if(rts=="drs"){add.constraint(lp.dea,c(rep(1,n)),indices=c(1:n),"<=",1)}
      
      # Input constraints
      for(i in 1:m){
        if(orientation=="i" && ncv[1,i]==0){add.constraint(lp.dea,c(xdata[,i],-xdata[k,i],1),indices=c(1:n,n+1,n+1+i),"=",0)}
        else if(orientation=="i" && ncv[1,i]==1){add.constraint(lp.dea,c(xdata[,i]-xdata[k,i],1),indices=c(1:n,n+1+i),"=",0)}
        else{add.constraint(lp.dea,c(xdata[,i],1),indices=c(1:n,n+1+i),"=",xdata[k,i])}
      }
      
      # Output constraints
      for(r in 1:s){
        if(orientation=="o" && ncv[1,m+r]==0){add.constraint(lp.dea,c(ydata[,r],-1*ydata[k,r],-1),indices=c(1:n,n+1,n+1+m+r),"=",0)}
        else if(orientation=="o" && ncv[1,m+r]==1){add.constraint(lp.dea,c(ydata[,r]-ydata[k,r],-1),indices=c(1:n,n+1+m+r),"=",0)}
        else{add.constraint(lp.dea,c(ydata[,r],-1),indices=c(1:n,n+1+m+r),"=",ydata[k,r])}
      }
      
      # External NDF
      if(!is.null(env)){for(j in 1:n){if(env[j,1]<env[k,1]){add.constraint(lp.dea,c(1),indices=c(j),"=",0)}}}
      
      # PPS for Super
      if(se==1){add.constraint(lp.dea,c(1),indices=c(k),"=",0)}
      
      # Bounds
      set.bounds(lp.dea,lower=c(rep(0,n),-Inf,rep(0,m+s)))  
      
      # Solve
      solve.lpExtPtr(lp.dea)
      
      # Get results
      results.efficiency[k]<-abs(get.objective(lp.dea))
      
      # Get results
      temp.p<-get.variables(lp.dea)
      results.lambda[k,]<-temp.p[1:n]
      results.xslack[k,]<-temp.p[(n+2):(n+1+m)]
      results.yslack[k,]<-temp.p[(n+1+m+1):(n+1+m+s)]
      temp.d<-get.dual.solution(lp.dea)
      results.vweight[k,]<-abs(temp.d[3:(2+m)])
      results.uweight[k,]<-abs(temp.d[(3+m):(2+m+s)])
      results.w[k,]<-temp.d[2]
      
      # Stage II
      if(exists("sg")){
        # Link previous solutions
        add.constraint(lp.dea,c(1),indices=c(n+1),"=",results.efficiency[k])
        
        # date sum
        if(sg=="max"){set.objfn(lp.dea,c(-date[1:n]),indices=c(1:n))}
        if(sg=="min"){set.objfn(lp.dea,c(date[1:n]),indices=c(1:n))}
        
        # slack sum max
        if(sg=="ssm"){set.objfn(lp.dea,c(rep(-1,m+s)),indices=c((n+2):(n+1+m+s)))}
        
        # solve
        solve.lpExtPtr(lp.dea)
        
        # get results
        temp.s<-get.variables(lp.dea)
        results.lambda[k,]<-temp.s[1:n]
        results.xslack[k,]<-temp.s[(n+2):(n+1+m)]
        results.yslack[k,]<-temp.s[(n+1+m+1):(n+1+m+s)]
      }
    }
    list(eff=results.efficiency,lambda=results.lambda,xslack=results.xslack,yslack=results.yslack,v=results.vweight,u=results.uweight,w=results.w)
  }
  
  # Loop for eff_r & eff_t
  for(i in 1:r){
    # Subset indices for each unique year
    if(i==1){s<-1}else{s<-till(d,unique(d)[i-1])+1}
    e<-till(d,unique(d)[i])
    x_t<-matrix(x[1:e,],nrow=e)
    y_t<-matrix(y[1:e,],nrow=e)
    d_t<-matrix(d[1:e,],nrow=e)
    
    # Run DEA
    if(i==r){temp<-dm.dea.internal(x_t,y_t,rts,orientation,0,sg,d_t,ncv,env,1,e)}
    else{temp<-dm.dea.internal(x_t,y_t,rts,orientation,0,sg,d_t,ncv,env,s,e)}
    
    # Save eff_r & eff_t
    if(i==r){eff_r[s:e,]<-temp$eff[s:e,];eff_t[1:e,]<-temp$eff[1:e,];lambda[1:e,1:e]<-temp$lambda[1:e,1:e]}
    else{eff_r[s:e,]<-temp$eff[s:e,]}
  }
  
  # Effective date
  if(ftype=="d"){
    for(i in 1:e){ed[i,1]<-sum(d[1:e,]*lambda[i,1:e]);sl[i,1]<-sum(lambda[i,1:e])}
    ed<-ed/sl
  }
  if(ftype=="s"){ed[,1]<-t}

  # RoC
  for(i in 1:e){
    if(round(eff_r[i,1],8)==1 && round(eff_t[i,1],8)!=1 && ed[i,1]>d[i,1]){
      if(orientation=="i"){roc[i,1]<-(1/eff_t[i,1])^(1/(ed[i,1]-d[i,1]))}
      if(orientation=="o"){roc[i,1]<-(eff_t[i,1])^(1/(ed[i,1]-d[i,1]))}
    }
  }
  
  # RoC filter
  roc[!is.na(roc[,1]) & roc[,1]>10,1]<-NA
  avgroc<-mean(roc,na.rm=TRUE)
  
  # RoC segmentation
  g_b<-array(0,c(n,1))
  g_b[!is.na(roc[,1]),1]<-1
  for(i in 1:e){
    # if(abs(eff_t[i,1]-1)<10^-9){local_roc[i,1]<-avgroc}
    if(sum(lambda[,i]*roc[,1],na.rm=TRUE)>0){local_roc[i,1]<-sum(lambda[,i]*roc[,1],na.rm=TRUE)/sum(lambda[,i]*g_b[,1],na.rm=TRUE)}
  }
  
  # Sort results back to original order
  eff_r<-matrix(c(eff_r[order(o),]),ncol=1)
  eff_t<-matrix(c(eff_t[order(o),]),ncol=1)
  lambda<-matrix(c(lambda[order(o),]),ncol=n)
  ed<-matrix(c(ed[order(o),]),ncol=1)
  roc<-matrix(c(roc[order(o),]),ncol=1)
  local_roc<-matrix(c(local_roc[order(o),]),ncol=1)
  
  results<-list(eff_r=eff_r,eff_t=eff_t,lambda_t=lambda,eft_date=ed,roc_past=roc,roc_local=local_roc,roc_avg=avgroc)
  return(results)
}
