#Objective function which we need to minimize 
#for the optimization problem
obj<- function(lam,u,ubar,Ti,rho,...){
  lam.t.u<- apply(u,1,crossprod,lam)
  #plot(lam.t.u)
  lam.t.ubar<- crossprod(lam,ubar)
  
  -mean(Ti*rho(lam.t.u, ...),na.rm = TRUE) + lam.t.ubar
}

#The first and Second derivative of the above objective function
derv.obj<- function(lam,Ti,rho1,u,ubar,...){
  lam.t.u <- apply(u,1,crossprod,lam)
  lam.t.ubar<- crossprod(lam,ubar)
  
  #get rho1(lam^Tu)*R/pi
  v<- numeric(length(Ti))
  v[Ti==1]<- rho1(lam.t.u, ...)[Ti==1]
  #get final answer
  -apply(apply(u,2,"*",v),2,mean,na.rm = TRUE)+ubar
}

derv2.obj<- function(lam,Ti,rho2,u,...){
  lam.t.u<- apply(u,1,crossprod,lam)
  
  #get rho(lam^Tu)
  v<- numeric(length(Ti))
  v[Ti==1]<- rho2(lam.t.u, ...)[Ti==1]
  
  #Get matrices for hessian
  mats<- sapply(1:nrow(u), function(i) tcrossprod(u[i,]),simplify = "array")
  mats2<- apply(mats,c(1,2),"*",v)
  #get final answer
  -apply(mats2,c(2,3),mean,na.rm = TRUE)
}

#The special case of CR family of functions

cr.rho<-function (v, theta) 
{
  v<-as.vector(v)
  if (theta == 0) {
    ans <- -exp(-v)
  }
  else if (theta == -1){
    ans <- suppressWarnings(log(1+v))
    ind <- is.na(ans)
    ans[ind] <- -Inf
  }
  else  {
    a <- -(1 - theta * v)^(1 + 1/theta)
    ans <- a/(theta + 1)
    ind <- is.na(ans)
    ans[ind] <- -Inf
  }
  return(ans)
}


#The first and second derivatives of the CR family functions

d.cr.rho<-  function (v, theta) 
{
  v<-as.vector(v)
  if (theta == 0) {
    ans <- exp(-v)
  }
  else if (theta == -1){
    ans <- 1/(1+v)
    ind <- is.na(ans)
    ans[ind] <- -Inf
  }
  else {
    ans <- (1 - theta * v)^(1/theta)
    ind <- is.na(ans)
    ans[ind] <- -Inf
  }
  return(ans)
}

dd.cr.rho<-   function (v, theta) 
{
  v<-as.vector(v)
  if (theta == 0) {
    a <- -exp(-v)
  }
  else if (theta == -1){
    a <- -1/(1+v)^2
    ind <- is.na(a)
    a[ind] <- -Inf
  }
  else {
    a <- -(1 - theta * v)^(1/theta - 1)
    ind <- is.na(a)
    a[ind] <- -Inf
  }
  a
}

###################################################################################
#The backtracking function for the main newton Raphson
backtrack<- function(alpha,beta,x.val,del.x,nabla.f,obj,u, ubar,Ti,rho, ...){
  u<- u
  step<- 1
  
  l.t.u<- apply(u,1,crossprod,x.val)
  f.x<- obj(x.val,u,ubar, Ti, rho,...)
  
  df.t.dx<- crossprod(nabla.f, del.x)
  #print(f.x)
  #print(df.t.dx)
  
  while(obj(x.val+step*del.x, u,ubar,Ti,rho,...) > f.x + alpha*step*df.t.dx ){
    step<- beta*step
  }
  return(step)
  
}

###################################################################################
