dm.ddf <-
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
  results.beta<-matrix(rep(NA,n*m),nrow=n,ncol=m) 
  results.gamma<-matrix(rep(NA,n*s),nrow=n,ncol=s) 
  results.xslack<-matrix(rep(NA,n*m),nrow=n,ncol=m) 
  results.yslack<-matrix(rep(NA,n*s),nrow=n,ncol=s)
  
  for (k in 1:n){   
    # Declare LP
    lp.ddf<-make.lp(0,n+n+m+s+m+s) # lambda+mu+beta+gamma+xslack+yslack
    
    # Set objective
    set.objfn(lp.ddf,c(rep(-1,m+s)),indices=c((n+n+1):(n+n+m+s)))
    
    # RTS
    if(rts=="vrs"){add.constraint(lp.ddf,c(rep(1,n*2)),indices=c(1:(n*2)),"=",1)}
    if(rts=="crs"){set.constr.type(lp.ddf,0,1)}
    if(rts=="irs"){add.constraint(lp.ddf,c(rep(1,n*2)),indices=c(1:(n*2)),">=",1)}
    if(rts=="drs"){add.constraint(lp.ddf,c(rep(1,n*2)),indices=c(1:(n*2)),"<=",1)}
    
    # Mu
    if(rts=="crs"||rts=="drs"||sum(w)==0){add.constraint(lp.ddf,c(rep(1,n)),indices=c((n+1):(n+n)),"=",0)}
    
    # Input constraints
    for(i in 1:m){add.constraint(lp.ddf,c(xdata[,i],xdata[,i],g[k,i],1),indices=c(1:(n+n),n+n+i,n+n+m+s+i),"=",xdata[k,i])}
    
    # Output constraints
    for(r in 1:s){
      if(w[1,r]==1){
        add.constraint(lp.ddf,c(ydata[,r],g[k,m+r]),indices=c(1:n,n+n+m+r),"=",ydata[k,r])
        add.constraint(lp.ddf,c(1),indices=c(n+n+m+s+m+r),"=",0)
      }else{add.constraint(lp.ddf,c(ydata[,r],-g[k,m+r],-1),indices=c(1:n,n+n+m+r,n+n+m+s+m+r),"=",ydata[k,r])}
      if(se==1){add.constraint(lp.ddf,c(ydata[,r],-1),indices=c(1:n,n+n+m+s+m+r),">=",0)}
    }
    
    # PPS for Super
    if(se==1){add.constraint(lp.ddf,c(1),indices=c(k),"=",0)}
    
    # Bounds
    set.bounds(lp.ddf,lower=c(rep(0,n+n+m+s+m+s)))  
    
    # Solve
    solve.lpExtPtr(lp.ddf)
    
    # Get results
    results.efficiency[k]<--1*get.objective(lp.ddf)
    
    # Get results
    temp.p<-get.variables(lp.ddf)
    results.lambda[k,]<-temp.p[1:n]
    results.mu[k,]<-temp.p[(n+1):(n+n)]
    results.beta[k,]<-temp.p[(n+n+1):(n+n+m)]
    results.gamma[k,]<-temp.p[(n+n+m+1):(n+n+m+s)]
    results.xslack[k,]<-temp.p[(n+n+m+s+1):(n+n+m+s+m)]
    results.yslack[k,]<-temp.p[(n+n+m+s+m+1):(n+n+m+s+m+s)]
    
    # Stage II
    if(exists("sg")){
      # Link previous solutions
      add.constraint(lp.ddf,c(rep(1,m+s)),indices=c((n+n+1):(n+n+m+s)),"=",results.efficiency[k])
      
      # date sum
      if(sg=="max"){set.objfn(lp.ddf,c(-date[1:n],-date[1:n]),indices=c(1:(n+n)))}
      if(sg=="min"){set.objfn(lp.ddf,c(date[1:n],date[1:n]),indices=c(1:(n+n)))}
      
      # slack sum max
      if(sg=="ssm"){set.objfn(lp.ddf,c(rep(-1,m+s)),indices=c((n+n+m+s+1):(n+n+m+s+m+s)))}
      
      # solve
      solve.lpExtPtr(lp.ddf)
      
      # get results
      temp.s<-get.variables(lp.ddf)
      results.lambda[k,]<-temp.s[1:n]
      results.mu[k,]<-temp.s[(n+1):(n+n)]
      results.beta[k,]<-temp.s[(n+n+1):(n+n+m)]
      results.gamma[k,]<-temp.s[(n+n+m+1):(n+n+m+s)]
      results.xslack[k,]<-temp.s[(n+n+m+s+1):(n+n+m+s+m)]
      results.yslack[k,]<-temp.s[(n+n+m+s+m+1):(n+n+m+s+m+s)]
    }
  }
  results<-list(eff=results.efficiency,lambda=results.lambda,mu=results.mu,beta=results.beta,gamma=results.gamma,xslack=results.xslack,yslack=results.yslack)
  return(results)
}
