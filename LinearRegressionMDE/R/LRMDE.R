#' Performs minimum distance estimation in linear regression model: Y=Xb + error
#'@param Y - Response variable in linear regression model
#'@param X - Explanatory variable in linear regression model
#'@return Returns betahat - Minimum distance estimator of b 
#'@examples
#'X <- matrix(c(1,1,3,4, 4,2), nrow=3, ncol=2)
#'Y <- c(1,-5, 8)
#'bhat <- LRMDE(Y,X) 
#'@references
#'[1] Koul, H. L (1985). Minimum distance estimation in linear regression with unknown error distributions. Statist. Probab. Lett., 3 1-8.
#'@references
#'[2] Koul, H. L (1986). Minimum distance estimation and goodness-of-fit tests in first-order autoregression. Ann. Statist., 14 1194-1213. 
#'@references
#'[3] Koul, H. L (2002). Weighted empirical process in nonlinear dynamic models. Springer, Berlin, Vol. 166
#'@seealso ARMDE
#'@export
#'@importFrom "stats" "optim"
#'@importFrom "utils" "read.csv"

LRMDE <- function(Y, X){
  
  DimMat <- dim(X)
  LengY <- length(Y)  
  if (is.null(DimMat) == TRUE ){
    message("X should be a matrix")
    stop		
  }else{
    nXRow <- DimMat[1]
    nXCol <- DimMat[2]
    
    if (nXCol != LengY){
      message("Dimension of X does not match dimension of Y")
      stop
    }
  }
  
  BetaOLS <- solve(t(X)%*%X)%*% (t(X)%*%Y)
  Tmin <- optim(BetaOLS, TLRLoss(Y, X))
  return(Tmin$par)
}


TLRLoss <- function(Y, X){

  DimMat <- dim(X)
  nXRow <- DimMat[1]
  nXCol <- DimMat[2]

  filepath = paste(system.file(package = "LinearRegressionMDE"),"/LegendreTable.csv", sep="")
  tbl <- read.csv(filepath, header = FALSE)
  x <- tbl[,1]
  w <- tbl[,2]
  nNode <- length(x)
  Dual <- function(t){
    fval = 0
    for (i in 1:nNode){
      
      for(k in 1:nXCol){
        fval <- fval +  (  ULRk_yt(k, x[i], t, Y, X) )^2 * w[i]
        
      }		
    }
    return(fval)
    
  }
  
  return(Dual)
}


ULRk_yt <- function(k, y, t, Y, X){
  
  DimMat <- dim(X)
  nXRow <- DimMat[1]
  nXCol <- DimMat[2]	
  
  
  fval <- 0;
  A <- (t(X)%*%X)^(-1/2);
  
  for(i in 1:nXRow){
    
    xi <- t(X[i,]) 
    ak <- A[,k]
    
    dik <- xi %*% ak
    
    fval <- fval + dik * (  ( Y[i]- xi %*% t  <= y ) -  ( -Y[i] + xi %*% t  < y ) ) 
  }
  fval <- fval/sqrt(nXRow)
  return(fval)
  
}
