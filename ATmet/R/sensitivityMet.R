sensitivityMet=function(model,x,y,nboot, method, conf)
{

  

  
  if (missing(nboot)){
    nboot <- 0
    warning("Confidence intervals for the sensitivity indices cannot be evaluated because no bootstrap replicates are used. Specify a number of bootstrap replicates in order to compute the confidence intervals",call.=TRUE)
  }
  if(missing(conf)){
    conf<-0.95
  }
  
  if(method!=("SRRC") && method!=("Sobol")){
    stop("The method for the calculation of the sensitivity indices should be either SRRC or Sobol")
  }
  
  
  if (method=="SRRC"){
    
    if (missing(model)&& missing(y)){
      stop("Please specify either a model function or the output quantity sample")
    }
    
    if (missing(y)){
      y=model(x)
    }
    
    if (length(y)!=dim(x)[1]) {
      stop("The length of the output quantity should be the same as the number of rows of the input quantities")  
    }
    SRRC <- src(x,y,rank=TRUE,nboot=nboot,conf)
    method <- "SRRC"
    S <- list( model=model, x=x, y=y, method=method, S1=SRRC$SRRC)
    return(S)
  }
  
  
  if (method=="Sobol"){
    #Separate the input sample into two independant samples
    if (dim(x)[1]%/%2==dim(x)[1]/2){
      x1 <- x[1:(dim(x)[1]/2),]
      x2 <- x[(dim(x)[1]/2+1):dim(x)[1],]}
    else {
      x1 <- x[1:(dim(x)[1]%/%2),]
      x2 <- x[(dim(x)[1]%/%2+1):(dim(x)[1]-1),]}
    
    Sob=sobol2007(model,x1,x2,nboot=nboot, conf)
    S=list(model=model, method=method, SI=Sob)
    return(S)
  }
}
