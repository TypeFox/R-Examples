
## Ipsen distance
##----------------------------------------
## ipsen <- function(object,...) UseMethod("ipsen")
ipsen <- function(object, ga=NULL, ...){
  if (is.null(ga)){
    if (object$tag == "undir"){
      optgamma <- optimal_gamma(object$N)
    } else {
      optgamma <- optimal_gamma_dir(object$N)
    }
  } else {
    optgamma <- ga
  }
  
  laplist <- object$L
  
  ## Check for parallelization
  n.cores <- NULL
  if (!is.na(match("n.cores",names(list(...)))))
    n.cores <- list(...)[["n.cores"]]

  verbose <- FALSE
  if (!is.na(match("verbose", names(list(...)))))
    verbose <- list(...)[["verbose"]]
  
  ## Should I use multiple cores or not?
  if(detectCores() >= 2 && (is.null(n.cores) || n.cores>1)){
    if (is.null(n.cores) || n.cores >= detectCores()){
      if (length(laplist) < detectCores()){
        n.cores <- length(laplist)
      } else {
        n.cores <- detectCores() - 1
      }
    }
    cl <- makeCluster(n.cores)
    
    ## Eval needed function on nodes
    clusterEvalQ(cl,{K <- function(mygamma,given_omega){
      return(1/integrate(lorentz,lower=0,upper=Inf,mygamma=mygamma,given_omega=given_omega, stop.on.error = FALSE)$value)
    }
                     rho <- function(omega, mygamma, ll){
                       ll[[2]]*lorentz(omega,mygamma,ll[[1]])
                     }
                     lorentz <- function(omega,mygamma,given_omega){
                       l <-0
                       for(i in 2:length(given_omega)){
                         l = l + mygamma/( (omega-given_omega[i])**2+mygamma**2)                          }
                       return(l)
                     }
                   })
    

    
    if (verbose)
      cat("Start computing eigenvalues with multiple cores\n")
    ## Actual computation of eigen-values/vectors
    ll <- clusterApply(cl,laplist,function(x,mygamma=optgamma,...){
      myomega <- sqrt(abs(round(spec(x),5)))
      myk <- K(mygamma,myomega)
      return(list(myomega,myk))
    })

    stopCluster(cl)
    
  } else {
    ## Computation on 1 CPU
    ll <- lapply(1:length(laplist),function(x,mygamma,laplist, ...){
      if (verbose)
        cat("Done",x,"/",length(laplist),"\n")
      aa <- laplist[[x]]
      myomega <- sqrt(abs(round(spec(aa),5)))
      myk <- K(mygamma,myomega)
      return(list(myomega,myk))
    }, mygamma=optgamma, laplist=laplist, ...)
  }
  mydistfun <- function(a,b, optgamma){
    integrand <- function(omega, mygamma, given_omega_G, given_omega_H){
      (rho(omega, optgamma,a)-rho(omega,optgamma,b))**2
    }
    tmp <- sqrt(integrate(integrand,lower=0,upper=Inf,mygamma=optgamma,given_omega_G=a[[1]],given_omega_H=b[[1]], stop.on.error=FALSE,rel.tol=.Machine$double.eps,subdivisions=1e4)$value)
    return(tmp)
  }
  if (verbose)
    cat("Start computing mutual distances\n")
  if (length(laplist) == 2){
    ## Compute distance between 2 adjacency matrices
    dist <- mydistfun(ll[[1]], ll[[2]], optgamma=optgamma)
    names(dist) <- "IM"
  } else {
    ## Compute mutual distances between all the matrices in the list
    idx <- combn(length(ll),2)
    tmpdist <- sapply(1:dim(idx)[2], function(x,ll,optgamma, idx, ...){
      if (verbose)
        cat("D(",idx[1,x],",", idx[2,x],")\n")
      mydistfun(ll[[idx[1,x]]], ll[[idx[2,x]]], optgamma)
    }, ll=ll, optgamma=optgamma, idx=idx)
    dist <- matrix(NA,ncol=length(ll), nrow=length(ll))
    dist[t(idx)] <- dist[t(idx)[,c(2,1)]] <- tmpdist
    diag(dist) <- 0
  }
  return(dist)
}


## Useful function for computing Ipsen distance
##--------------------------------------------------
spec <- function(mm){
  sort(eigen(mm)$values)
}

## D2 - ipsen02evolutionary
lorentz <- function(omega,mygamma,given_omega){
  l <-0
  for(i in 2:length(given_omega)){
    l = l + mygamma/( (omega-given_omega[i])**2+mygamma**2)                          }
  return(l)
}

K <- function(mygamma,given_omega){
  return(1/integrate(lorentz,lower=0,upper=Inf,mygamma=mygamma,given_omega=given_omega, stop.on.error = FALSE)$value)
}

rho <- function(omega, mygamma, ll){
  ll[[2]]*lorentz(omega,mygamma,ll[[1]])
}

ipsen_minus_one <- function(g,n){
  return(sqrt(
    1/(pi*g) +
      1/(2*g*(atan(sqrt(n)/g)+pi/2)**2)*(pi/2+ (sqrt(n)/g)/(1+(sqrt(n)/g)**2)+atan(sqrt(n)/g))
    -4*(pi-(g/sqrt(n))*log(1/(1+(sqrt(n)/g)**2))+atan(sqrt(n)/g))/
      (pi*g*(4+(sqrt(n)/g)**2)*(atan(sqrt(n)/g)+pi/2))
  )
         -1)
}

optimal_gamma <- function(n){
  return(uniroot(ipsen_minus_one,c(0.01,1),n=n,maxiter=100000,tol=.Machine$double.eps)$root)
}


##--------------------------------------------------
## DIRECTED GRAPH DISTANCE
##--------------------------------------------------

ipsen_minus_one_dir  <- function(g,n){
  return(ZZ(g)^2*MM(0,g)+WW(n,g)^2*MM(n-2,g)+WW(n,g)^2*MM(n,g)+WWp(n,g)^2*MM(2*n-2,g)
         -2*ZZ(g)*WW(n,g)*LL(0,n-2,g)-2*ZZ(g)*WW(n,g)*LL(0,n,g)-2*ZZ(g)*WWp(n,g)*LL(0,2*n-2,g)
         +2*WW(n,g)*WW(n,g)*LL(n-2,n,g) +2*WW(n,g)*WWp(n,g)*LL(n-2,2*n-2,g) +2*WW(n,g)*WWp(n,g)*LL(n,2*n-2,g)-1)
}

optimal_gamma_dir <- function(n){
  return(uniroot(ipsen_minus_one_dir,c(0.01,1),n=n,maxiter=100000,tol=.Machine$double.eps)$root)
}


MM <- function(T,g){
  return(1/2*(g^2*atan(sqrt(T)/g) + T*atan(sqrt(T)/g) + sqrt(T)*g)/(g^5 +T*g^3) + 1/4*pi/g^3)
}

LL  <- function(T,U,g){
  return(-log(g^2 + U)/((4*g^2 + T + 3*U)*sqrt(T) - (4*g^2 + 3*T + U)*sqrt(U)) + log(g^2 + T)/((4*g^2 + T + 3*U)*sqrt(T) - (4*g^2 + 3*T + U)*sqrt(U)) + pi/(4*g^3 + T*g - 2*sqrt(T)*sqrt(U)*g + U*g) + atan(sqrt(T)/g)/(4*g^3 + T*g - 2*sqrt(T)*sqrt(U)*g + U*g) + atan(sqrt(U)/g)/(4*g^3 + T*g - 2*sqrt(T)*sqrt(U)*g + U*g))
}

ZZ  <- function(g){
  return(2*g/pi)
}

WW  <- function(N,g){
  WWW  <- (2*N-1)*(pi/2)+(N-1)*(atan(sqrt(N-2)/g) +atan(sqrt(N)/g)) + atan(sqrt(2*N-2)/g)
  return(g*(N-1)/WWW)
}

WWp  <- function(N,g){
  return(WW(N,g)/(N-1))
}


