
hessian <- function (func, x, method="Richardson",
                              method.args=list(), ...) UseMethod("hessian")

hessian.default <- function(func, x, method="Richardson",
      method.args=list(), ...){  

 if(1!=length(func(x, ...)))
       stop("Richardson method for hessian assumes a scalar valued function.")

 if(method=="complex"){ # Complex step hessian
   args <- list(eps=1e-4, d=0.1, 
      zero.tol=sqrt(.Machine$double.eps/7e-7), r=4, v=2)
   args[names(method.args)] <- method.args
   # the CSD part of this uses eps=.Machine$double.eps
   # but the jacobian is Richardson and uses method.args
   return(jacobian(func=function(fn, x, ...){grad(func=fn, x=x, 
           method="complex", method.args=list(eps=.Machine$double.eps), ...)}, 
         x=x, fn=func, method.args=args, ...))
   } 
 else if(method != "Richardson")  stop("method not implemented.")

   args <- list(eps=1e-4, d=0.1, zero.tol=sqrt(.Machine$double.eps/7e-7), 
                r=4, v=2, show.details=FALSE) # default
   args[names(method.args)] <- method.args
   D <- genD(func, x, method=method, method.args=args, ...)$D
   if(1!=nrow(D)) stop("BUG! should not get here.")
   H <- diag(NA,length(x))
   u <- length(x)
   for(i in 1:length(x))
      for(j in 1:i){ 
         u <- u + 1
         H[i,j] <- D[,u]
         }

   H <- H + t(H)
   diag(H) <- diag(H)/2
   H
  }


#######################################################################

#               Bates & Watts   D matrix calculation

#######################################################################

genD <- function(func, x, method="Richardson",
                   method.args=list(), ...)UseMethod("genD")

genD.default <- function(func, x, method="Richardson",
      method.args=list(), ...){
  #   additional cleanup by Paul Gilbert (March, 2006)
  #   modified substantially by Paul Gilbert (May, 1992)
  #    from original code by Xingqiao Liu,   May, 1991.
  
  #  This function is not optimized for S speed, but is organized in 
  # the same way it could be (was) implemented in C, to facilitate checking.
  
  #  v  reduction factor for Richardson iterations. This could
  #	 be a parameter but the way the formula is coded it is assumed to be 2.
 
    if(method != "Richardson")  stop("method not implemented.")
    args <- list(eps=1e-4, d=0.0001, zero.tol=sqrt(.Machine$double.eps/7e-7),
              r=4, v=2) # default
    args[names(method.args)] <- method.args
    d <- args$d
    r <- args$r
    v <- args$v
    if (v!=2) stop("The current code assumes v is 2 (the default).")	 

    f0 <- func(x, ...)   #f0 is the value of the function at x.

    n <- length(x)  #  number of parameters (theta)
    h0 <- abs(d*x) + args$eps * (abs(x) < args$zero.tol)
    D <- matrix(0, length(f0),(n*(n + 3))/2)
    #length(f0) is the dim of the sample space
    #(n*(n + 3))/2 is the number of columns of matrix D.( first
    #	der. & lower triangle of Hessian)
    Daprox <-  matrix(0,length(f0),r) 
    Hdiag  <-  matrix(0,length(f0),n)
    Haprox <-  matrix(0,length(f0),r)

    for(i in 1:n){    # each parameter  - first deriv. & hessian diagonal
          h <-h0
          for(k in 1:r){  # successively reduce h 
             f1 <- func(x+(i==(1:n))*h, ...)
             f2 <- func(x-(i==(1:n))*h, ...) 
             #f1 <- do.call("func",append(list(x+(i==(1:n))*h), func.args))
             #f2 <- do.call("func",append(list(x-(i==(1:n))*h), func.args))
             Daprox[,k] <- (f1 - f2)  / (2*h[i])    # F'(i) 
             Haprox[,k] <- (f1-2*f0+f2)/ h[i]^2     # F''(i,i) hessian diagonal
             h <- h/v	  # Reduced h by 1/v.
             }
          for(m in 1:(r - 1))
             for ( k in 1:(r-m)){
                Daprox[,k]<-(Daprox[,k+1]*(4^m)-Daprox[,k])/(4^m-1)
        	Haprox[,k]<-(Haprox[,k+1]*(4^m)-Haprox[,k])/(4^m-1)
                }
          D[,i] <- Daprox[,1]
          Hdiag[,i] <- Haprox[,1]
          }	 
    u <- n
  
    for(i in 1:n){   # 2nd derivative  - do lower half of hessian only
       for(j in 1:i){
          u <- u + 1
             if (i==j) D[,u] <- Hdiag[,i]
             else {
               h <-h0
               for(k in 1:r){  # successively reduce h 
        	  f1 <- func(x+(i==(1:n))*h + (j==(1:n))*h, ...)
        	  f2 <- func(x-(i==(1:n))*h - (j==(1:n))*h, ...)
        	  Daprox[,k]<- (f1 - 2*f0 + f2 -
        			   Hdiag[,i]*h[i]^2 - 
        			   Hdiag[,j]*h[j]^2)/(2*h[i]*h[j])  # F''(i,j)  
        	  h <- h/v     # Reduced h by 1/v.
        	  }
               for(m in 1:(r - 1))
        	  for ( k in 1:(r-m))
        	     Daprox[,k]<-(Daprox[,k+1]*(4^m)-Daprox[,k])/(4^m-1)
               D[,u] <- Daprox[,1]
               }
          }  
       }
  D <- list(D=D, p=length(x), f0=f0, func=func, x=x, d=d,
            method=method, method.args=args)# Darray constructor (genD.default)
  class(D) <- "Darray"
  invisible(D)
  }
  
