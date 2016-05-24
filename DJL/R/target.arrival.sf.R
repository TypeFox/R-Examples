target.arrival.sf <-
function(xdata,ydata,date,t,rts,g,w=NULL,sg="ssm",ftype="d"){
  
  # Initial checks
  if(is.na(match(rts,c("crs","vrs","irs","drs")))){stop('rts must be "crs", "vrs", "irs", or "drs".')}
  if(is.na(match(sg,c("ssm","max","min")))){stop('sg must be "ssm", "max", or "min".')}
  if(is.na(match(ftype,c("d","s")))){stop('ftype must be either "d" or "s".')}
  if(t<=min(date)){stop('t is earlier than dataset.')}
  if(max(date)<=t){stop('t is later than dataset.')}
  
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
  eff_t<-array(NA,c(n,1))
  eff_t_gm<-array(NA,c(n,1)) # Geometric mean for equi-ratio
  lambda<-array(NA,c(n,n))
  ed<-array(NA,c(n,1))
  sl<-array(NA,c(n,1))
  roc_ind<-array(NA,c(n,1))
  arrival_avg<-array(NA,c(n,1))
  arrival_seg<-array(NA,c(n,1))
  
  # Subset index
  till<-function(x,y){
    t<-0
    while(x[t+1]<=y&&t<nrow(x)){t<-t+1}
    return(t)
  }

  #Subset till t
  e<-till(d,t)
  x_t<-matrix(x[1:e,],nrow=e)
  y_t<-matrix(y[1:e,],nrow=e)
  g_t<-matrix(g[1:e,],nrow=e)
  
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
  
  # Loop for eff_t
  for(i in (e+1):n){
    # Subset till t + each target
    x_f<-rbind(x_t,x[i,])
    y_f<-rbind(y_t,y[i,])
    g_f<-rbind(g_t,g[i,])
    # Run SF
    temp<-dm.sf.internal(x_f,y_f,rts,g_f,w,se=1,sg,d,e+1,e+1)
    # Save eff_t
    eff_t[i,]<-temp$eff[e+1,]
    lambda[i,1:e]<-temp$lambda[e+1,1:e]
  }

  # Effective date
  if(ftype=="d"){for(i in (e+1):n){ed[i,1]<-sum(d[1:e,]*lambda[i,1:e])/sum(lambda[i,1:e])}}
  if(ftype=="s"){ed[,1]<-t}
  
  # Calc iRoC
  roc<-roc.sf(xdata,ydata,date,t,rts,g,w,sg,ftype)
  roc_local<-roc$roc_local;roc_local_bi<-ifelse(is.na(roc_local),0,1);roc_avg<-roc$roc_avg
  for(i in (e+1):n){roc_ind[i,1]<-sum(roc_local[1:e,]*lambda[i,1:e],na.rm=TRUE)/sum(lambda[i,1:e]*roc_local_bi[1:e,])}

  # Arrival target
  for(i in (e+1):n){
    if(abs(eff_t[i,1])>10^-9 && abs(eff_t[i,1])!=Inf && eff_t[i,1]<0){
      eff_t_gm[i,1]<-((1-eff_t[i,1])/(1+eff_t[i,1]))^0.5 # Geometric mean for equi-ratio
      arrival_avg[i,1]<-ed[i,1]+log(eff_t_gm[i,1],exp(1))/log(roc_avg,exp(1))
      arrival_seg[i,1]<-ed[i,1]+log(eff_t_gm[i,1],exp(1))/log(ifelse(roc_ind[i,1]>0,roc_ind[i,1],roc_avg),exp(1))
      #arrival_seg[i,1]<-ed[i,1]+log(1/eff_t[i,1],exp(1))/log(roc_ind[i,1],exp(1))
    }
  }
  results<-list(eff_t=eff_t,lambda_t=lambda,eft_date=ed,roc_avg=roc_avg,roc_local=roc_local,roc_ind=roc_ind,arrival_avg=arrival_avg,arrival_seg=arrival_seg)
  return(results)    
}
