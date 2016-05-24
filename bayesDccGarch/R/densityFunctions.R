
dssnorm <-
function(x, gamma=rep(1,length(x)), log=FALSE){
  if(any(gamma <= 0)) stop("The gamma parameters must be greater than zero.")

  value = .C("dssged_R", x=as.double(x), 
    gamma=as.double(gamma), delta=as.double(2.0), 
    k=as.integer(length(x)), islog = as.integer(log), 
    value=numeric(1))$value
  
  return(value)
}


dsst <-
function(x, gamma=rep(1,length(x)), nu=10, log=FALSE){
  if(nu < 2) stop("The nu parameter must be greater than 2.")
  if(any(gamma <= 0)) stop("The gamma parameters must be greater than zero.")

  value = .C("dsst_R", x=as.double(x),
    gamma=as.double(gamma), v=as.double(nu), 
    k=as.integer(length(x)), islog = as.integer(log), 
    value=numeric(1))$value
  
  return(value)
}


dssged <-
function(x, gamma=rep(1,length(x)), delta=2, log=FALSE){
  if(delta <= 0) stop("The delta parameter must be greater than zero.")
  if(any(gamma <= 0)) stop("The gamma parameters must be greater than zero.")
  
  value = .C("dssged_R", x=as.double(x), 
    gamma=as.double(gamma), delta=as.double(delta), 
    k=as.integer(length(x)), islog = as.integer(log), 
    value=numeric(1))$value

  return(value)
}

