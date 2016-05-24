## (C) Jeffrey S. Racine July 22 2011

## glp.model.matrix is a modified version of the polym() function
## (stats) combined with the tensor.prod.model.matrix function in
## mgcv. The function accepts a vector of degrees and provides a
## generalized polynomial with varying polynomial order. This can be
## more parsimonious than the tensor product model matrix commonly
## found in the spline literature while retaining solid approximation
## capabilities and will be better conditioned than the tensor product
## (see commented illustration below).

## X is a list of model matrices, from which a generalized local
## polynomial model matrix is to be produced (first column ones).

## Note - if this is used for modeling and derivatives are required,
## one would pass in the derivative model matrix for the variable(s)
## whose derivative is required and pass matrices of zeros of
## identical dimension to the originals for all other variables.

glp.model.matrix <- function(X) {

  k <-length(X)
  dimen.list <- list()
  dimen <- numeric()
  for(i in 1:k) {
    dimen[i] <- ncol(X[[i]])
    dimen.list[[i]] <- 0:dimen[i]
  }
  z <- do.call("expand.grid", dimen.list, k)
  s <- rowSums(z)
  ind <- (s > 0) & (s <= max(dimen))
  z <- z[ind, ,drop=FALSE]
  if(!all(dimen==max(dimen))) {
    for(j in 1:length(dimen)) {
      d <- dimen[j]
      if((d < max(dimen)) & (d > 0)) {
        s <- rowSums(z)
        d <- (s > d) & (z[,j,drop=FALSE]==matrix(d,nrow(z),1,byrow=TRUE))
        z <- z[!d, ]
      }
    }
  }

  res <- cbind(1, X[[1]])[, 1 + z[, 1]]
  if(k > 1) for (i in 2:k) res <- res * cbind(1, X[[i]])[, 1 + z[, i]]
  return(matrix(res,nrow=NROW(X[[1]])))

}

##> set.seed(42)
##> n <- 1000
##> 
##> degree <- 10
##> nbreak <- 2
##> 
##> X1 <- gsl.bs(runif(n),degree=degree,nbreak=nbreak)
##> X2 <- gsl.bs(runif(n),degree=degree,nbreak=nbreak)
##> X3 <- gsl.bs(runif(n),degree=degree,nbreak=nbreak)
##> 
##> X <- list()
##> X[[1]] <- X1
##> X[[2]] <- X2
##> X[[3]] <- X3
##> 
##>   k<-length(X)
##>   dimen.list <- list()
##>   for(i in 1:k) dimen.list[[i]] <- 0:ncol(X[[i]])
##> 
##> B.glp <- glp.model.matrix(X)
##> B.tp <- tensor.prod.model.matrix(X)
##> 
##> dim(B.glp)
##[1] 1000  285
##> dim(B.tp)
##[1] 1000 1000
##> 
##> all.equal(matrix(B.glp),matrix(B.tp))
##[1] "Attributes: < Component 1: Mean relative difference: 2.508772 >"
##[2] "Numeric: lengths (285000, 1000000) differ"                      
##> 
##> rcond(t(B.glp)%*%B.glp)
##[1] 8.338926e-13
##> rcond(t(B.tp)%*%B.tp)
##[1] 5.92409e-21

