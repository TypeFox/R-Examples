  # grad case 1 and 2 are special cases of jacobian, with a scalar rather than
  #   vector valued function. Case 3 differs only because of the interpretation
  #   that the vector result is a scalar function applied to each argument, and the
  #   thus the result has the same length as the argument.
  #   The code of grad could be consolidated to use jacobian.
  #   There is also some duplication in genD.

############################################################################

#    functions for gradient calculation

############################################################################

grad <- function (func, x, method="Richardson", side=NULL, 
   method.args=list(), ...) UseMethod("grad")

grad.default <- function(func, x, method="Richardson", side=NULL,
      method.args=list(), ...){
  # modified by Paul Gilbert from code by Xingqiao Liu.
  # case 1/ scalar arg, scalar result (case 2/ or 3/ code should work)
  # case 2/ vector arg, scalar result (same as special case jacobian)
  # case 3/ vector arg, vector result (of same length, really 1/ applied multiple times))
  f <- func(x, ...)
  n <- length(x)	 #number of variables in argument

  if (is.null(side)) side <- rep(NA, n)
  else {
       if(n != length(side)) 
          stop("Non-NULL argument 'side' should have the same length as x")
       if(any(1 != abs(side[!is.na(side)]))) 
          stop("Non-NULL argument 'side' should have values NA, +1, or -1.")
       }

  case1or3 <- n == length(f)

  if((1 != length(f)) & !case1or3)
  	 stop("grad assumes a scalar valued function.")

  if(method=="simple"){
    #  very simple numerical approximation
    args <- list(eps=1e-4) # default
    args[names(method.args)] <- method.args

    side[is.na(side)] <- 1
    eps <- rep(args$eps, n) * side

    if(case1or3) return((func(x+eps, ...)-f)/eps) 

    # now case 2
    df <- rep(NA,n)
    for (i in 1:n) {
      dx <- x
      dx[i] <- dx[i] + eps[i] 
      df[i] <- (func(dx, ...) - f)/eps[i]
     }
    return(df)
    } 
  else if(method=="complex"){ # Complex step gradient
    if (any(!is.na(side))) stop("method 'complex' does not support non-NULL argument 'side'.")
    eps <- .Machine$double.eps
    v <- try(func(x + eps * 1i, ...))
    if(inherits(v, "try-error")) 
      stop("function does not accept complex argument as required by method 'complex'.")
    if(!is.complex(v)) 
      stop("function does not return a complex value as required by method 'complex'.")
   
    if(case1or3) return(Im(v)/eps) 
    # now case 2
    h0 <- rep(0, n)
    g  <- rep(NA, n)
    for (i in 1:n) {
      h0[i] <- eps * 1i
      g[i] <- Im(func(x+h0, ...))/eps 
      h0[i]  <- 0
      }
    return(g)
    } 
  else if(method=="Richardson"){
    args <- list(eps=1e-4, d=0.0001, zero.tol=sqrt(.Machine$double.eps/7e-7), r=4, v=2, show.details=FALSE) # default
    args[names(method.args)] <- method.args
    d <- args$d
    r <- args$r
    v <- args$v
    show.details <- args$show.details
    a <- matrix(NA, r, n) 
    #b <- matrix(NA, (r - 1), n)
  
    #  first order derivatives are stored in the matrix a[k,i], 
    #  where the indexing variables k for rows(1 to r), i for columns (1 to n),
    #  r is the number of iterations, and n is the number of variables.
  

    h <- abs(d*x) + args$eps * (abs(x) < args$zero.tol)
    pna <- (side == 1)  & !is.na(side) # double these on plus side
    mna <- (side == -1) & !is.na(side) # double these on minus side

    for(k in 1:r)  { # successively reduce h		    
       ph <- mh <- h
       ph[pna] <- 2 * ph[pna] 
       ph[mna] <- 0           
       mh[mna] <- 2 * mh[mna] 
       mh[pna] <- 0           

       if(case1or3)  a[k,] <- (func(x + ph, ...) -  func(x - mh, ...))/(2*h)
       else for(i in 1:n)  {
    	 if((k != 1) && (abs(a[(k-1),i]) < 1e-20)) a[k,i] <- 0 #some func are unstable near zero
    	 else  a[k,i] <- (func(x + ph*(i==seq(n)), ...) - 
    			  func(x - mh*(i==seq(n)), ...))/(2*h[i])
    	 }
       if (any(is.na(a[k,]))) stop("function returns NA at ", h," distance from x.")
       h <- h/v     # Reduced h by 1/v.
       }	

   if(show.details)  {
        cat("\n","first order approximations", "\n")		
        print(a, 12)
    }

  #------------------------------------------------------------------------
  # 1 Applying Richardson Extrapolation to improve the accuracy of 
  #   the first and second order derivatives. The algorithm as follows:
  #
  #   --  For each column of the derivative matrix a,
  #	  say, A1, A2, ..., Ar, by Richardson Extrapolation, to calculate a
  #	  new sequence of approximations B1, B2, ..., Br used the formula
  #
  #	     B(i) =( A(i+1)*4^m - A(i) ) / (4^m - 1) ,  i=1,2,...,r-m
  #
  #		N.B. This formula assumes v=2.
  #
  #   -- Initially m is taken as 1  and then the process is repeated 
  #	 restarting with the latest improved values and increasing the 
  #	 value of m by one each until m equals r-1
  #
  # 2 Display the improved derivatives for each
  #   m from 1 to r-1 if the argument show.details=T.
  #
  # 3 Return the final improved  derivative vector.
  #-------------------------------------------------------------------------
  
    for(m in 1:(r - 1)) {	  
       a <- (a[2:(r+1-m),,drop=FALSE]*(4^m)-a[1:(r-m),,drop=FALSE])/(4^m-1)
       if(show.details & m!=(r-1) )  {
  	  cat("\n","Richarson improvement group No. ", m, "\n") 	  
  	  print(a[1:(r-m),,drop=FALSE], 12)
  	}
     }
  return(c(a))
  } else stop("indicated method ", method, "not supported.")
}
  

jacobian <- function (func, x, method="Richardson", side=NULL,
                              method.args=list(), ...) UseMethod("jacobian")

jacobian.default <- function(func, x, method="Richardson", side=NULL,
      method.args=list(), ...){
  f <- func(x, ...)
  n <- length(x)	 #number of variables.

  if (is.null(side)) side <- rep(NA, n)
  else {
       if(n != length(side)) 
          stop("Non-NULL argument 'side' should have the same length as x")
       if(any(1 != abs(side[!is.na(side)]))) 
          stop("Non-NULL argument 'side' should have values NA, +1, or -1.")
       }

  if(method=="simple"){
    #  very simple numerical approximation
    args <- list(eps=1e-4) # default
    args[names(method.args)] <- method.args

    side[is.na(side)] <- 1
    eps <- rep(args$eps, n) * side

    df <-matrix(NA, length(f), n)
    for (i in 1:n) {
      dx <- x
      dx[i] <- dx[i] + eps[i] 
      df[,i] <- (func(dx, ...) - f)/eps[i]
     }
    return(df)
    } 
  else if(method=="complex"){ # Complex step gradient
    if (any(!is.na(side))) stop("method 'complex' does not support non-NULL argument 'side'.")
    # Complex step Jacobian
    eps <- .Machine$double.eps
    h0  <-  rep(0, n)
    h0[1] <- eps * 1i
    v <- try(func(x+h0, ...))
    if(inherits(v, "try-error")) 
      stop("function does not accept complex argument as required by method 'complex'.")
    if(!is.complex(v)) 
      stop("function does not return a complex value as required by method 'complex'.")
  
    h0[1]  <- 0
    jac <- matrix(NA, length(v), n)
    jac[, 1] <- Im(v)/eps
    if (n == 1) return(jac)
    for (i in 2:n) {
      h0[i] <- eps * 1i
      jac[, i] <- Im(func(x+h0, ...))/eps 
      h0[i]  <- 0
      }
    return(jac)
    } 
  else if(method=="Richardson"){
    args <- list(eps=1e-4, d=0.0001, zero.tol=sqrt(.Machine$double.eps/7e-7), 
                r=4, v=2, show.details=FALSE) # default
    args[names(method.args)] <- method.args
    d <- args$d
    r <- args$r
    v <- args$v 	  
    a <- array(NA, c(length(f),r, n) )
  
    h <- abs(d*x) + args$eps * (abs(x) < args$zero.tol)
    pna <- (side == 1)  & !is.na(side) # double these on plus side
    mna <- (side == -1) & !is.na(side) # double these on minus side

    for(k in 1:r)  { # successively reduce h		 
       ph <- mh <- h
       ph[pna] <- 2 * ph[pna] 
       ph[mna] <- 0           
       mh[mna] <- 2 * mh[mna] 
       mh[pna] <- 0           

       for(i in 1:n)  {
    	 a[,k,i] <- (func(x + ph*(i==seq(n)), ...) -  
     		     func(x - mh*(i==seq(n)), ...))/(2*h[i])
    	 #if((k != 1)) a[,(abs(a[,(k-1),i]) < 1e-20)] <- 0 #some func are unstable near zero
    	 }
       h <- h/v     # Reduced h by 1/v.
       }     

   for(m in 1:(r - 1)) {	  
       a <- (a[,2:(r+1-m),,drop=FALSE]*(4^m)-a[,1:(r-m),,drop=FALSE])/(4^m-1)
     }
  # drop second dim of a, which is now 1 (but not other dim's even if they are 1
  return(array(a, dim(a)[c(1,3)]))  
  } else stop("indicated method ", method, "not supported.")
}

