dm.dea <-
function(xdata,ydata,rts,orientation,se=0,sg="ssm",date=NULL,ncv=NULL,env=NULL){

  # Initial checks
  if(is.na(match(rts,c("crs","vrs","irs","drs")))){stop('rts must be "crs", "vrs", "irs", or "drs".')}
  if(is.na(match(orientation,c("i","o")))){stop('orientation must be either "i" or "o".')}
  if(is.na(match(se,c(0,1)))){stop('se must be either 0 or 1.')}
  if(is.na(match(sg,c("ssm","max","min")))){stop('sg must be "ssm", "max", or "min".')}
  
  # Load library
  # library(lpSolveAPI)  

  # Parameters
  xdata<-as.matrix(xdata);ydata<-as.matrix(ydata);if(!is.null(date))date<-as.matrix(date) # format input data as matrix
  n<-nrow(xdata); m<-ncol(xdata); s<-ncol(ydata)
  if(is.null(ncv)) ncv<-matrix(c(0),ncol=m+s) else ncv<-as.matrix(ncv)
  if(!is.null(env))env<-as.matrix(env)
  
  # Data frames
  results.efficiency<-matrix(rep(NA,n),nrow=n,ncol=1)
  results.lambda<-matrix(rep(NA,n^2),nrow=n,ncol=n)
  results.xslack<-matrix(rep(NA,n*m),nrow=n,ncol=m) 
  results.yslack<-matrix(rep(NA,n*s),nrow=n,ncol=s) 
  results.vweight<-matrix(rep(NA,n*m),nrow=n,ncol=m) 
  results.uweight<-matrix(rep(NA,n*s),nrow=n,ncol=s) 
  results.w<-matrix(rep(NA,n),nrow=n,ncol=1)

  for (k in 1:n){
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
  results<-list(eff=results.efficiency,lambda=results.lambda,xslack=results.xslack,yslack=results.yslack,v=results.vweight,u=results.uweight,w=results.w)
  return(results)
}
