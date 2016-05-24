dm.sf <-
function(xdata,ydata,rts,g,w=NULL,se=0,sg="ssm",date=NULL){

  # Initial checks
  if(is.na(match(rts,c("crs","vrs","irs","drs")))){stop('rts must be "crs", "vrs", "irs", or "drs".')}
  if(is.na(match(se,c(0,1)))){stop('se must be either 0 or 1.')}
  if(is.na(match(sg,c("ssm","max","min")))){stop('sg must be "ssm", "max", or "min".')}
  
  # Load library
  # library(lpSolveAPI)
  
  # Parameters
  xdata<-as.matrix(xdata);ydata<-as.matrix(ydata);g<-as.matrix(g);if(!is.null(date))date<-as.matrix(date) # format input data as matrix
  n<-nrow(xdata); m<-ncol(xdata); s<-ncol(ydata)
  if(is.null(w)) w<-matrix(c(0),ncol=s) else w<-as.matrix(w)
  
  # Data frames
  results.efficiency<-matrix(rep(NA,n),nrow=n,ncol=1)
  results.lambda<-matrix(rep(NA,n^2),nrow=n,ncol=n)
  results.mu<-matrix(rep(NA,n^2),nrow=n,ncol=n)
  results.xslack<-matrix(rep(NA,n*m),nrow=n,ncol=m) 
  results.yslack<-matrix(rep(NA,n*s),nrow=n,ncol=s) 
  
  for (k in 1:n){   
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
  results<-list(eff=results.efficiency,lambda=results.lambda,mu=results.mu,xslack=results.xslack,yslack=results.yslack)
  return(results)
}
