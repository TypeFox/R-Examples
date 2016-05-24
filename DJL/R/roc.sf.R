roc.sf <-
function(xdata,ydata,date,t,rts,g,w=NULL,sg="ssm",ftype="d"){
  
  # Initial checks
  if(is.na(match(rts,c("crs","vrs","irs","drs")))){stop('rts must be "crs", "vrs", "irs", or "drs".')}
  if(is.na(match(sg,c("ssm","max","min")))){stop('sg must be "ssm", "max", or "min".')}
  if(is.na(match(ftype,c("d","s")))){stop('ftype must be either "d" or "s".')}
  if(t<=min(date)){stop('t is earlier than dataset.')}
  if(max(date)<t){stop('t is later than dataset.')}
  
  # Parameters
  xdata<-as.matrix(xdata);ydata<-as.matrix(ydata);date<-as.matrix(date);g<-as.matrix(g) # format input data as matrix
  n<-nrow(xdata); m<-ncol(xdata); s<-ncol(ydata)
  if(is.null(w)) w<-matrix(c(0),ncol=s) else w<-as.matrix(w)
  o<-matrix(c(1:n),ncol=1) # original data order
  
  # Sort data ascending order
  x<-matrix(c(xdata[order(date),]),ncol=m)
  y<-matrix(c(ydata[order(date),]),ncol=s)
  d<-matrix(c(date[order(date),]),ncol=1)
  g<-matrix(c(g[order(date),]),ncol=m+s)
  o<-matrix(c(o[order(date),]),ncol=1)
  
  # Data frames
  eff_r<-array(NA,c(n,1))
  eff_t<-array(NA,c(n,1))
  eff_t_gm<-array(NA,c(n,1)) # Geometric mean for equi-ratio
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
  
  # SF internal function 
  dm.sf.internal<-function(xdata,ydata,rts,g,w,se=0,sg,date,a,z){
    
    # Load library
    # library(lpSolveAPI)
    
    # Parameters
    n<-nrow(xdata); m<-ncol(xdata); s<-ncol(ydata)
    if(is.null(w)){w<-matrix(c(0),ncol=s)}
    
    # Data frames
    results.efficiency<-matrix(rep(NA,n),nrow=n,ncol=1)
    results.lambda<-matrix(rep(NA,n^2),nrow=n,ncol=n)
    results.mu<-matrix(rep(NA,n^2),nrow=n,ncol=n)
    results.xslack<-matrix(rep(NA,n*m),nrow=n,ncol=m) 
    results.yslack<-matrix(rep(NA,n*s),nrow=n,ncol=s) 
    
    for (k in a:z){   
      # Declare LP
      lp.sf<-make.lp(0,n+n+1+m+s) # lambda+mu+efficiency+xslack+yslack
      
      # Set objective
      set.objfn(lp.sf,c(-1),indices=c(n+n+1))
      
      # RTS
      if(rts=="vrs"){add.constraint(lp.sf,c(rep(1,n*2)),indices=c(1:(n*2)),"=",1)}
      if(rts=="crs"){set.constr.type(lp.sf,0,1)}
      if(rts=="irs"){add.constraint(lp.sf,c(rep(1,n*2)),indices=c(1:(n*2)),">=",1)}
      if(rts=="drs"){add.constraint(lp.sf,c(rep(1,n*2)),indices=c(1:(n*2)),"<=",1)}
      
      # Mu
      if(rts=="crs"||rts=="drs"||sum(w)==0){add.constraint(lp.sf,c(rep(1,n)),indices=c((n+1):(n+n)),"=",0)}
      
      # Input constraints
      for(i in 1:m){add.constraint(lp.sf,c(xdata[,i],xdata[,i],g[k,i],1),indices=c(1:n,(n+1):(n+n),n+n+1,n+n+1+i),"=",xdata[k,i])}
      
      # Output constraints
      for(r in 1:s){
        if(w[1,r]==1){
          add.constraint(lp.sf,c(ydata[,r],ydata[,r],g[k,m+r]),indices=c(1:n,(n+1):(n+n),n+n+1),"=",ydata[k,r])
          add.constraint(lp.sf,c(1),indices=c(n+n+1+m+r),"=",0)
        }else{add.constraint(lp.sf,c(ydata[,r],-g[k,m+r],-1),indices=c(1:n,n+n+1,n+n+1+m+r),"=",ydata[k,r])}
        if(se==1){add.constraint(lp.sf,c(ydata[,r],-1),indices=c(1:n,n+n+1+m+r),">=",0)}
      }
      
      # PPS for Super
      if(se==1){add.constraint(lp.sf,c(1,1),indices=c(k,n+k),"=",0)}
      
      # Bounds
      set.bounds(lp.sf,lower=c(rep(0,n+n),-Inf,rep(0,m+s)))  
      
      # Solve
      solve.lpExtPtr(lp.sf)
      
      # Get results
      results.efficiency[k]<--1*get.objective(lp.sf)

      # Get results
      temp.p<-get.variables(lp.sf)
      results.lambda[k,]<-temp.p[1:n]
      results.mu[k,]<-temp.p[(n+1):(n+n)]
      results.xslack[k,]<-temp.p[(n+n+2):(n+n+1+m)]
      results.yslack[k,]<-temp.p[(n+n+1+m+1):(n+n+1+m+s)]
      
      # Stage II
      if(exists("sg")){
        # Link previous solutions
        add.constraint(lp.sf,c(1),indices=c(n+n+1),"=",results.efficiency[k])
        
        # date sum
        if(sg=="max"){set.objfn(lp.sf,c(-date[1:n],-date[1:n]),indices=c(1:n,(n+1):(n+n)))}
        if(sg=="min"){set.objfn(lp.sf,c(date[1:n],date[1:n]),indices=c(1:n,(n+1):(n+n)))}
        
        # slack sum max
        if(sg=="ssm"){set.objfn(lp.sf,c(rep(-1,m+s)),indices=c((n+n+2):(n+n+1+m+s)))}
        
        # solve
        solve.lpExtPtr(lp.sf)
        
        # get results
        temp.s<-get.variables(lp.sf)
        results.lambda[k,]<-temp.s[1:n]
        results.mu[k,]<-temp.s[(n+1):(n+n)]
        results.xslack[k,]<-temp.s[(n+n+2):(n+n+1+m)]
        results.yslack[k,]<-temp.s[(n+n+1+m+1):(n+n+1+m+s)]
      }
    }
    list(eff=results.efficiency,lambda=results.lambda,mu=results.mu,xslack=results.xslack,yslack=results.yslack)
  }

  # Loop for eff_r & eff_t
  for(i in 1:r){
    # Subset indices for each unique year
    if(i==1){s<-1}else{s<-till(d,unique(d)[i-1])+1}
    e<-till(d,unique(d)[i])
    x_t<-matrix(x[1:e,],nrow=e)
    y_t<-matrix(y[1:e,],nrow=e)
    d_t<-matrix(d[1:e,],nrow=e)
    g_t<-matrix(g[1:e,],nrow=e)
    
    # Run SF
    if(i==r){temp<-dm.sf.internal(x_t,y_t,rts,g_t,w,0,sg,d_t,1,e)}
    else{temp<-dm.sf.internal(x_t,y_t,rts,g_t,w,0,sg,d_t,s,e)}
    
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
    if(round(eff_r[i,1],8)==0 && round(eff_t[i,1],8)!=0 && ed[i,1]>d[i,1]){
      eff_t_gm[i,1]<-((1+eff_t[i,1])/(1-eff_t[i,1]))^0.5 # Geometric mean for equi-ratio
      roc[i,1]<-(eff_t_gm[i,1])^(1/(ed[i,1]-d[i,1]))
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
