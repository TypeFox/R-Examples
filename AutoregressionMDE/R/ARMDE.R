#' Performs minimum distance estimation in autoregressive model
#'@param X : vector of n observed value
#'@param AR_Order : oder of the autoregressive model
#'@return returns minimum distance estimators of the parameter in the autoregressive model
#'@examples
#'X <- rnorm(10, mean=0, sd=1)
#'AR_Order <- 2
#'rhohat<-ARMDE(X,AR_Order)
#'@references
#'[1] Koul, H. L (1985). Minimum distance estimation in linear regression with unknown error distributions. Statist. Probab. Lett., 3 1-8.
#'@references
#'[2] Koul, H. L (1986). Minimum distance estimation and goodness-of-fit tests in first-order autoregression. Ann. Statist., 14 1194-1213. 
#'@references
#'[3] Koul, H. L (2002). Weighted empirical process in nonlinear dynamic models. Springer, Berlin, Vol. 166
#'@importFrom "stats" "optim"
#'@importFrom "utils" "read.csv"
#'@export
#'@seealso LRMDE

ARMDE <- function(X, AR_Order){
  
  nLength <- length(X)
  
  if(nLength<=AR_Order){
    message("Length of vector X should be greater than AR_Order.")
    stop
  }
  
  Xres <- rep(0, times=(nLength-AR_Order))
  tempvec <- rep(0, times= AR_Order*(nLength-AR_Order) ) 
  Xexp <- matrix( tempvec, nrow = (nLength-AR_Order), ncol = AR_Order)
  
  for (i in 1:(nLength-AR_Order) ) {
    Xres[i] <- X[nLength - (i-1)]
    for (j in 1:AR_Order){
      Xexp[i,j] <- X[nLength-(i+j-1) ]
        
    }
  }               
  
  tempdet = det(  t(Xexp) %*% Xexp )
  if (  tempdet < 0.01 ){
    rho0 <- 0.5*rep(1, times = AR_Order)  
  }else{
    rho0 <- solve(t(Xexp)%*%Xexp)%*% (t(Xexp)%*%Xres)
  }
    
  Tmin <- optim(rho0, TLoss(X, AR_Order))
  return(Tmin$par)
}


TLoss <- function(X, AR_Order){
  
  filepath <- paste(system.file(package = "AutoregressionMDE"),"/LegendreTable.csv", sep="")
  tbl <- read.csv(filepath, header = FALSE)
  x <- tbl[,1]
  w <- tbl[,2]
  nNode <- length(x)
  
  Dual <- function(t){
    fval <- 0
    for (i in 1:nNode){
      
      for(k in 1:AR_Order){
        fval <- fval +  (  Uk_yt(k, x[i], t, X, AR_Order) )^2 * w[i]
        
      }		
    }
    return(fval)
    
  }
  
  return(Dual)
  
  
}

Uk_yt <- function(k, y, t, X, AR_Order){
  
  nLength <- length(X)
  fval <- 0
  
  for(i in 1:nLength){
    if(i-k<=0){
      Xik <- 0
    }else{
      Xik <- X[i-k]
    }
    tXi1 <- 0
    
    for(j in 1:AR_Order){
      if(i-j<=0){
        Xij <- 0
      }else{
        Xij <- X[i-j]
      }	
      tXi1 <- tXi1 + t[j]*Xij;	
    }
    fval <- fval + Xik*(  (X[i]-tXi1<=y  ) -  ( -X[i] + tXi1 < y  )  )
  }
  
  fval <- fval/sqrt(nLength)
  return(fval)
}