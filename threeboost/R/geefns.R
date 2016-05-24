#' GEE estimating functions
#' 
#' Internal functions for computing the GEE. Should generally not be called by user.
#' @export
#' @param Y Vector of (correlated) outcomes
#' @param X Matrix of predictors
#' @param b Vector of coefficients
#' @param mu.Y Mean function
#' @param g.Y Link function (inverse of mean function)
#' @param v.Y Variance function
#' @param aux Auxiliary function for coputing (co)variance parameters
#' @param id Vector of cluster IDs
#' @param uid Vector of unique subject IDs 
#' @param rows.indivs List of rows of \code{X} corresponding to each subject ID
#' @param corstr Working correlation structure
ee.GEE <- function(Y,X,b,mu.Y,g.Y,v.Y,aux,id=1:length(Y),uid=sort(unique(id)),
                   rows.indivs=lapply(uid,function(j) { which(id==j)}),corstr="ind") {

  aux.par <- aux(Y=Y,X=X,b=b,id=id)
  a <- aux.par[1]
  phi <- aux.par[2]
  b <- c(aux.par[3],b)
  X <- cbind(rep(1,nrow(X)),X)
  
  eta <- t(X%*%b)
  
  if(corstr=="exch") {
    ##print("Setting up exchangeable working correlation matrix")    
    indivMats <- lapply(rows.indivs,function(rs) { 
      clust.size <- length(rs)
      sqA.i <- diag(sqrt(v.Y(eta[rs])))
      R <- stats::toeplitz(c(1,rep(a,clust.size-1)))
      phi*solve(sqA.i%*%R%*%sqA.i)
    })    
    Vinv <- as.matrix(Matrix::bdiag(indivMats))
  } else {
    Vinv <- phi*diag(as.vector(1/v.Y(eta)))
  }
  
  contrib.indiv <- function(ind,b) {
    eta.indiv <- X[ind,]%*%b
    if(length(ind)==1) {
      A <- v.Y(eta.indiv)
      cont <- A*X[ind,]*Vinv[ind,ind]*(Y[ind] - mu.Y(eta.indiv))
    }
    else {
      A <- diag(as.vector(v.Y(eta.indiv)))
      cont <- t(A %*% X[ind,]) %*% Vinv[ind,ind] %*% (Y[ind] - mu.Y(eta.indiv))
    }
    cont
  }
  
  L <- lapply(rows.indivs,contrib.indiv,b=b)
  
  contribs <- rowSums(do.call("cbind",L))
  return(contribs[-1]) ## Don't return the element of the intercept
}

#' @describeIn ee.GEE
#' @export
ee.GEE.aux <- function(Y,X,b,mu.Y,g.Y,v.Y,id=1:length(Y),uid=sort(unique(id)),
                       rows.indivs=lapply(uid,function(j) { which(id==j)})) {
  
  eta <- t(X%*%b)
  pearson.resids <- (Y - mu.Y(eta)) / sqrt(v.Y(eta))
  
  p <- sum(b!=0) - 1 ## Really, this p should be number of covariates in model. We subtract 1 to exclude intercept
  
  num <- sum(unlist(lapply(rows.indivs,function(rs) {
    combs <- combn(rs,2)
    sum(pearson.resids[combs[1,]]*pearson.resids[combs[2,]])    
  })))
  
  ## This function could be moved "outside" because it doesn't need to be recomputed every time.
  denom <- sum(unlist(lapply(rows.indivs,function(rs) {
    n.i <- length(rs)
    (n.i*(n.i-1))/2
  })))- p
  
  phi <- 1/(sum(pearson.resids^2)/(length(Y)-p))
  
  alpha <- phi*num/denom ## This gives an estimate of pairwise correlation alpha
  
  ## (Re)-estimate the intercept
  mu.hat <- mu.Y(eta)
  b0 <- mean(g.Y(mu.hat) - as.vector(X[,-1]%*%b[-1]))
  
  return(c(alpha,phi,b0))
  
}

expit <- function(x) { exp(x)/(1+exp(x))}

mu.Lin <- function(eta){eta}
mu.Bin <- function(eta){expit(eta)}
mu.Pois <- function(eta){exp(eta)}

g.Lin <- function(m){m}
g.Bin <- function(m){log(m/(1-m))}
g.Pois <- function(m){ log(m) }

v.Lin <- function(eta){rep(1,length(eta))}
v.Bin <- function(eta){ expit(eta)/(1+exp(eta))}
v.Pois <- function(eta) { exp(eta) }

ee.GEELin <- function(...) { ee.GEE(...,mu.Y=mu.Lin,g.Y=g.Lin,v.Y=v.Lin) }
ee.GEEBin <- function(...) { ee.GEE(...,mu.Y=mu.Bin,g.Y=g.Bin,v.Y=v.Bin )}
ee.GEEPois <- function(...) {ee.GEE(...,mu.Y=mu.Pois,g.Y=g.Pois,v.Y=v.Pois)}

ee.GEELin.aux <- function(...) { ee.GEE.aux(...,mu.Y=mu.Lin,g.Y=g.Lin,v.Y=v.Lin) }
ee.GEEBin.aux <- function(...) { ee.GEE.aux(...,mu.Y=mu.Bin,g.Y=g.Bin,v.Y=v.Bin) }
ee.GEEPois.aux <- function(...) { ee.GEE.aux(...,mu.Y=mu.Pois,g.Y=g.Pois,v.Y=v.Pois)}
