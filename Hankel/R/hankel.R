
.hankel.statisticExp <- function(x,y,m,n,lambda=1.0,standardized=FALSE)
{
  
  if(standardized==FALSE) out <- .C("teststatistic", as.double(x),as.double(y),as.integer(m),as.integer(n),as.double(lambda),result=as.double(1))
  
  else{
    
    out <- .C("teststatisticStandardized", as.double(x),as.double(y),as.integer(m),as.integer(n),as.double(lambda),as.double(mean(c(x,y))),result=as.double(1))
    
  }  
  
  return(out$result)
}


.hankel.critValExp <- function(x,y,m,n,lambda,standardized,bootstrapSampleSize,sigLevel)
{

  xy<-c(x,y)
  T_mnValues<-c()
  
  for(i in 1:bootstrapSampleSize)
  {
    xresampled<-sample(xy,m,replace=TRUE)
    yresampled<-sample(xy,n,replace=TRUE)
    T_mnValues[i]<-.hankel.statisticExp(xresampled,yresampled,m,n,lambda,standardized)
  }
 
  return(quantile(T_mnValues,probs=sigLevel,type=1))
  
}


.hankel.pValExp <- function(x,y,m,n,lambda,standardized,bootstrapSampleSize,statistic)
{
  
  xy<-c(x,y)
  T_mnValues<-c()
  
  for(i in 1:bootstrapSampleSize)
  {
    xresampled<-sample(xy,m,replace=TRUE)
    yresampled<-sample(xy,n,replace=TRUE)
    T_mnValues[i]<-.hankel.statisticExp(xresampled,yresampled,m,n,lambda,standardized)
  }

  return(mean(T_mnValues>=statistic))
  
}


hankel.test <- function(x,y,standardized=FALSE,replicates=500,calcCritVal=FALSE,calcPVal=TRUE,sigLevel=0.95,probMeasure="exp",params=list(lambda=1.0))
{
  
  if(!is.vector(x)  ||  !is.vector (y)){ stop("x,y must be a vector")}
  if(!is.numeric(x) || !is.numeric(y)){ stop("x,y must be numeric")}
  if(length(x)!=length(x[x>=0]) || length(y)!=length(y[y>=0]) )  {stop("x,y must be non-negative")}
  
  
  if(!is.numeric(sigLevel)){ stop("sigLevel must be numeric")}
  if(!(sigLevel>0) || !(sigLevel<1)){"sigLevel must be within (0,1)"}  
  if(!is.numeric(sigLevel)){ stop("sigLevel must be numeric")}
  
  if(probMeasure!="exp")  {stop("Unknown probMeasure set")} #To be changed when new probability measures are being implemented.
  
  if(probMeasure=="exp"){
    
    if(!is.numeric(params$lambda) ) { stop("lambda must be numeric")}
    if(!(params$lambda>0) ) { stop("lambda must be non-negative")}
    
    Summary<-list(
               samplesizes=c(length(x),length(y)),
               statistic="not computed",
               sigLevel=sigLevel,
               calcPVal=calcPVal,
               pValue="not computed",
               critValue="not computed",
               calcCritVal=calcCritVal,
               replicates=as.integer(replicates),
               standardized=standardized,
               params=params,
               probMeasure="exp"
              )
    
    Summary$statistic<-.hankel.statisticExp(x,y,Summary$samplesizes[1],Summary$samplesizes[2],params$lambda,Summary$standardized)
    
    if(calcCritVal) Summary$critValue<-.hankel.critValExp(x,y,Summary$samplesizes[1],Summary$samplesizes[2],params$lambda,Summary$standardized,Summary$replicates,Summary$sigLevel)
    if(calcPVal) Summary$pValue <- .hankel.pValExp(x,y,Summary$samplesizes[1],Summary$samplesizes[2],params$lambda,Summary$standardized,Summary$replicates,Summary$statistic)
  }  
    
  class(Summary) <- "hankeltest"
  return(Summary)
}


print.hankeltest<-function(x,...) {
  
  cat("X-samplesize:",x$samplesizes[1],"\n")
  cat("Y-samplesize:",x$samplesizes[2],"\n")
  cat("Probability measure:",x$probMeasure,", lambda =",x$params$lambda,"\n")
  cat("Standardized version:" ,x$standardized,"\n")
  cat("Statistic:",x$statistic,"\n")
  
  if (x$calcCritVal || x$calcPVal ){
    
      cat("Number of bootstrap replicates:",x$replicates,"\n")
    
      if (x$calcCritVal)  {
      
        cat("Significance level:",x$sigLevel,"\n")
        cat("Critical value for significance level:",x$critValue,"\n")
            
     }  
    
     if (x$calcPVal)  {
       
       cat("P-value:",x$pValue,"\n")
       
     }  
     
     
} 
  
  invisible(x)
}




