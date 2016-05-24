################################################################################
#####    The fista methods for MRSP                                        #####
################################################################################
#####    Author: Wolfgang Pößnecker                                        #####
#####    Last modified: 27.05.2014, 20:10                                  #####
################################################################################

##### functions for computing the l2norm
## always returns a vector, even if the input was a matrix or array
l2normraw <- cmpfun(function(u) sqrt(sum(c(u)^2)))
l2norm <- cmpfun(function(u){sqrt(length(c(u)) * sum(c(u)^2))})
l2norminv <- cmpfun(function(u){sqrt(sum(c(u)^2) / length(c(u)))})
l2norm.x <- cmpfun(function(u){if(!is.matrix(u)) u <- as.matrix(u)
 if(is.matrix(u)) sqrt((length(c(u))-ncol(u)) * sum(c(u)^2))
 else stop("error in l2norm: non-matrix argument")})
l2norminv.x <- cmpfun(function(u){if(!is.matrix(u)) u <- as.matrix(u)
 if(is.matrix(u)){sqrt(sum(c(u)^2) / (length(c(u))-ncol(u)))}else{
 stop("error in l2norminv: non-matrix argument")}})            
## a little helper function to replace apply in the computation of the penalty 
## for ordinary (unstructured) lasso penalty types.
rowL2norm <- cmpfun(function(x){d<-nrow(x);out<-numeric(d);for(i in 1:d){out[i]<-l2norm(x[i,])};out})
rowL2norminv <- cmpfun(function(x){d<-nrow(x);out<-numeric(d);for(i in 1:d){out[i]<-l2norminv(x[i,])};out})
rowL2normraw <- cmpfun(function(x){d<-nrow(x);out<-numeric(d);for(i in 1:d){out[i]<-l2normraw(x[i,])};out})


## a helper function that computes the projection onto the L2-norm ball:
l2project <- cmpfun(function(u, w, lambda, eps, groupweight){
 norm <- l2normraw(u)
 if(lambda*w*norm/eps*groupweight > 1){
  u <- u/norm
 }else{
  u <- lambda*w*groupweight/eps*u
 }                
 u
})

## a helper function that computes the projection onto the L_infinity-norm ball:
lmaxproject <- cmpfun(function(u){                                                                            ## hier dieselben argumente einbauen wie in l2project!
 u <- pmin(u, 1)
 u <- pmax(u, -1)
 u
})


## a helper function that computes the squared spectral norm (aka largest eigen-
## value) of the matrix w*R. for adjacent- or all-pairwise-FusedLasso, it can be
## computed analytically, but not if weights are used (ie w != rep(1, nrow(R)).
spectralnorm <- cmpfun(function(A, squared = T){
 if(!is.matrix(A)) stop("A must be a matrix in 'spectralnorm'")
 val <- eigen(crossprod(A), symmetric=F, only.values=T)$values[1]
 if(!squared) val <- sqrt(val)
 val
})



##### some helper functions for the case of symmetric side constraints
## a function to compute the sum that has to be zero after applying tresh
constrfun <- cmpfun(function(x, y, lambda){
u <- t(t(y)-x)
st <- 1-lambda/sqrt(diag(tcrossprod(u)))
z <- c(pmax(st,0))*u
return(abs(colSums(z)))})

## a function to adjust the output of tresh if the constraints are not satisfied
correctfun <- cmpfun(function(x, y, lambda, ind, tol=1e-3){
 testfun <- function(x){sum(constrfun(x,y,lambda))}
 for(i in 1:length(x)){
  if(ind[i]){
   testfuni <- function(u){
    x[i] <- u
    eval(call("testfun", x))
   }
   x[i] <- as.numeric(optimize(testfuni, c(x[i]-3, x[i]+3), tol = tol)[1])
  }
 }
 x    
}) 
 

symfun <- cmpfun(function(x, lambda){
 if(!is.matrix(x))
   stop("something went wrong with the structuring of the coef-objects!")
 if(ncol(x) == 1){
  m <- (min(x) + max(x)) / 2
  indlist <- abs(x - m) > lambda
  if(!any(indlist)){out <- m}else{
   signlist <- sign(m - x)
   out <- sum(x[indlist] + signlist[indlist]*lambda[indlist]) / length(x[indlist])
  }
  out
 }else{
  u <- t(x)
  m <- (-maxRow(-u) + maxRow(u)) / 2
  u <- t(u-m)
  indlist <- vector(length=nrow(u))
  for(i in 1:length(indlist)) indlist[i] <- l2normraw(u[i,]) > lambda[i]
  out <- numeric(ncol(x))
  if(all(indlist == F)){out <- m}else{
   for(i in 1:length(out)){
    signlist <- (m[i] - x[,i])/l2normraw((m[i] - x[,i]))
    out[i] <- sum(x[indlist,i] + signlist[indlist]*lambda[indlist]) / length(x[indlist,i])
   }
  }
  ind <- constrfun(out, x, lambda) > 0.1
  if(any(ind)) out <- correctfun(out, x, lambda, ind)
  #if(sum(constrfun(out, x, lambda)) > 0.1) out <- correctfun(out, x, lambda)  
  out
 }
 return(out)  
})



##### proximal operator for global predictors
## l2norm and generic soft thresholding function for models with ref category
treshrefx <- cmpfun(function(u, lambda, w, nn){
 if(nn) u <- pmax(u, 0)
 w <- sqrt(length(c(u)) - ncol(u)) * w                                                
 st <- 1 - lambda * w / l2normraw(u)                                            
 u <- max(st, 0) * u
 u
})

## soft thresholding function for ordinary lasso penalties.
treshlrefx <- cmpfun(function(u, lambda, w, nn){
 #w <- sqrt(ncol(u)) * w
 #st <- 1 - lambda * w / sqrt(diag(tcrossprod(u)))
 #u <- c(pmax(st, 0)) * u
 if(nn) u <- pmax(u, 0)
 dp <- ncol(u)
 dn <- nrow(u)
 norms <- rep(0, dn)
 out <- u
 u <- matrix(lassoC(as.double(u), as.double(lambda), as.double(w), as.integer(dp),
                as.integer(dn), as.double(norms), as.double(out))[[7]], nrow = dn)
 u
})
 
## sparse thresholding for sparse group lasso. input u must be a Q*P matrix.
sptreshrefx <- cmpfun(function(u, lambda, lambdasp, w, wsp, nn){
 u <- treshlrefx(u, lambdasp, wsp, nn)
 u <- treshrefx(u, lambda, w, nn)
 u
})

## P-spline tresholding.
Psplinetreshrefx <- cmpfun(function(u, lambda, Omega, tD){
 dp <- ncol(u)
 #dn <- nrow(u)
 w <- ncol(tD)
 #for(i in seq(dn)){
 # u[i,] <- solve(diag(dp) + lambda * w * Omega) %*% u[i,,drop=T]
 #}
 u <- u %*% solve(diag(dp) + lambda * w * Omega)
 u
})

## column-wise P-spline thresholding
colPsplinetreshrefx <- cmpfun(function(u, lambda, R){
 dp <- ncol(u)
 dn <- nrow(u)
 w <- dp
 RtR <- crossprod(R)
 u <- solve(diag(dn) + lambda * w * RtR) %*% u
 u
})



## l2norm and generic soft thresholding function for models with sym. constr.
## note for self: the weight (singular!) doesnt have to be used for centering!!!
treshsymx <- cmpfun(function(u, lambda, w, nn){
 w <- sqrt(length(c(u)) - ncol(u)) * w
 u <- t(t(u) - treshmean(u))
 if(nn) u <- pmax(u, 0)
 st <- 1 - lambda * w / l2normraw(u)
 u <- max(st, 0) * u
 u
})

## soft thresholding function for ordinary lasso penalties.
treshlsymx <- cmpfun(function(u, lambda, w, nn){
 #w <- sqrt(ncol(u)) * w
 ####u <- t(t(u) - symfun(u, w*lambda))
 #st <- 1 - lambda * w / sqrt(diag(tcrossprod(u)))
 #u <- c(pmax(st, 0)) * u
 if(nn) u <- pmax(u, 0)
 dp <- ncol(u)
 dn <- nrow(u)
 norms <- rep(0, dn)
 out <- u
 u <- matrix(lassoC(as.double(u), as.double(lambda), as.double(w), as.integer(dp),
                as.integer(dn), as.double(norms), as.double(out))[[7]], nrow = dn)
 u
})
 
## sparse thresholding for sparse group lasso. input u must be a Q*P matrix.
sptreshsymx <- cmpfun(function(u, lambda, lambdasp, w, wsp, nn){
 u <- treshlsymx(u, lambdasp, wsp, nn)
 u <- treshsymx(u, lambda, w, nn)
 u
})

## P-spline tresholding.
Psplinetreshsymx <- cmpfun(function(u, lambda, Omega, tD){
 dp <- ncol(u)
 #dn <- nrow(u)
 w <- ncol(tD)
 #for(i in seq(dn)){
 # u[i,] <- solve(diag(dp) + lambda * w * Omega) %*% u[i,,drop=T]
 #}
 u <- u %*% solve(diag(dp) + lambda * w * Omega)
 u
})

## column-wise P-spline thresholding
colPsplinetreshsymx <- cmpfun(function(u, lambda, R){
 dp <- ncol(u)
 dn <- nrow(u)
 w <- dp
 RtR <- crossprod(R)
 u <- solve(diag(dn) + lambda * w * RtR) %*% u
 u
})

###################################################################################################
## sparse thresholding for sparse group lasso combined with L1 fusion based on ADMM.
spFGLtreshsymx <- cmpfun(function(v, lambdal2, lambdal1, lambdaF, wl2, wl1, wF, R, nn, trace=F){

 ## algorithmic tuning parameter
 rho <- 0.1 #1/lambdaF #0.1
 rhoscale <- 1 + 1/rho

 ## initialize stuff
 ## adjust weights for problem size
 wl2a <- sqrt(length(c(v)) - ncol(v)) * wl2
 wl1a <- sqrt(ncol(v)) * wl1
 wFa <- sqrt(ncol(v)) * wF

 ## a function for computing the objective function of spFGLtreshnonex
 objective <- function(x){0.5*sum((x - v)^2) + lambdal2*wl2a*l2normraw(x) +
                          lambdal1*sum(wl1a*rowL2normraw(x)) + lambdaF*sum(wFa*rowL2normraw(R%*%x))}

 ## algorithmic tuning parameters
 rel.tol <- 5e-3
 iter.max <- 200

 ## precompute and store t(R) and R'R
 tR <- t(R)
 #RtR <- crossprod(R)

 ## store the correct tresholding function belonging to the sparse group lasso part of the penalty
 if(lambdal2 * lambdal1 > 0) xtresh <- function(x) sptreshsymx(x, lambdal2/rhoscale, lambdal1/rhoscale, wl2, wl1, nn)
 else if(lambdal2 == 0) xtresh <- function(x) treshlsymx(x, lambdal1/rhoscale, wl1, nn)
 else if(lambdal1 == 0) xtresh <- function(x) treshsymx(x, lambdal2/rhoscale, wl2, nn)

 ## store the correct tresholding function for the z-update
 ztresh <- function(z) treshlsymx(z, lambdaF, wF, nn)

 ## initialize the algorithmic quantities
 if(nn) v <- pmax(v, 0)
 x <- xtresh(v)
 z <- ztresh(R%*%x)
 u <- matrix(0, nrow=nrow(R), ncol=ncol(v))               # this is the scaled residual

 ## prepare for the loop
 iter.count <- 0
 d.fn <- 1e30
 d.fn.counter <- 0
 d.fn.converged <- F
 fn.val <- fn.val.old<- 1e30
 primal.converged <- dual.converged <- F

 #### the main ADMM loop ####
 while(!d.fn.converged || !primal.converged || !dual.converged){
  iter.count <- iter.count + 1
  if(iter.count > iter.max){
   if(trace){
    warning(paste("maximum number of iterations reached in spGFLtreshnonex"))
   }
   break
  }

  z.old <- z
  fn.val.old <- fn.val

  linpoint <- x  -  rho * tR %*% (R%*%x - z + u)
  xpoint <- (rho*v + linpoint)/(1 + rho)

  x <- xtresh(xpoint)                           # x-update
  Rx <- R%*%x

  zpoint <- Rx + u
  z <- ztresh(zpoint)                           # z-update

  primal.res <- Rx - z                          # primal residual
  u <- u + primal.res                           # u-update

  dual.res <- tR%*%(z - z.old)                  # dual residual                 #/ lambdaF ### hier eventuell durch lambdaF oder rho teilen?

  if(l2normraw(primal.res) <= rel.tol * max(c(l2normraw(Rx), l2normraw(z)))){
   primal.converged <- T
  }
  #if(l2normraw(dual.res) <= rel.tol * l2normraw(tR%*%u)){  ### hier eventuell noch rechte Seite mal lambdaF oder rho nehmen?
   dual.converged <- T
  #}

  if(!d.fn.converged){
   fn.val <- objective(x)                          # value of the objective
   d.fn <- abs(fn.val - fn.val.old)/(rel.tol/10 + abs(fn.val))
   if(d.fn < rel.tol) d.fn.converged <- T
  }
 } ## end while(!d.fn.converged ...)
 x
})

## l2norm and generic soft thresholding function for models with no constraint.
## (i.e. identifiability is obtained solely by the penalty.)
treshnonex <- cmpfun(function(u, lambda, w, nn){
 if(nn) u <- pmax(u, 0)
 w <- sqrt(length(c(u))) * w                                                
 st <- 1 - lambda * w / l2normraw(u)                                            
 u <- max(st, 0) * u
 u
})

## soft thresholding function for ordinary lasso penalties.
treshlnonex <- cmpfun(function(u, lambda, w, nn){
 #w <- sqrt(ncol(u)) * w
 #st <- 1 - lambda * w / sqrt(diag(tcrossprod(u)))
 #u <- c(pmax(st, 0)) * u
 if(nn) u <- pmax(u, 0)
 dp <- ncol(u)
 dn <- nrow(u)
 norms <- rep(0, dn)
 out <- u
 u <- matrix(lassoC(as.double(u), as.double(lambda), as.double(w), as.integer(dp),
                    as.integer(dn), as.double(norms), as.double(out))[[7]], nrow = dn)
 u
})
 
## sparse thresholding for sparse group lasso. input u must be a Q*P matrix.
sptreshnonex <- cmpfun(function(u, lambda, lambdasp, w, wsp, nn){
 u <- treshlnonex(u, lambdasp, wsp, nn)
 u <- treshnonex(u, lambda, w, nn)
 u
})

## P-spline tresholding.
Psplinetreshnonex <- cmpfun(function(u, lambda, Omega, tD){
 dp <- ncol(u)
 #dn <- nrow(u)
 w <- ncol(tD)
 #for(i in seq(dn)){
 # u[i,] <- solve(diag(dp) + lambda * w * Omega) %*% u[i,,drop=T]
 #}
 u <- u %*% solve(diag(dp) + lambda * w * Omega)
 u
})

## column-wise P-spline thresholding
colPsplinetreshnonex <- cmpfun(function(u, lambda, R){
 dp <- ncol(u)
 dn <- nrow(u)
 w <- dp
 RtR <- crossprod(R)
 u <- solve(diag(dn) + lambda * w * RtR) %*% u
 u
})


## group Bridge penalty:
proxfunGbridge <- cmpfun(function(u, l, w=matrix(1,nrow=nrow(u),ncol=ncol(u)), nn=F){
 w <- pmin(w, 1e38)                  # necessary for numerical reasons...
 u.orig <- u
 w.orig <- w
 u <- c(u)
 w <- c(w)
 f <- function(x){ 1/2*crossprod(c(x)-c(u.orig)) + 2*sqrt(l)*sqrt(sum(abs(c(w)*c(x)))) }   # sqrt(l*p) ?
 if(nn) u <- pmax(u, 0)
 w.ind <- which(w > 0)               # index of positive weights
 wt <- w[w.ind]
 p <- length(w.ind)
 if(p < 2) stop("the case 'length(which(w > 0)) < 2' is not implemented yet in proxfunGbridge")
 v <- abs(u)[w.ind]/wt
 myorder <- order(v, decreasing=T)   # sort in decreasing order
 v <- v[myorder]
 wt <- wt[myorder]
 vwsum <- cumsum(v*wt^2)
 wsum <- cumsum(wt^2)
 taus <- numeric(p)
 for(i in p:1){                      # for each interval [v[i+1], v[i]), find the possible roots via cubic equation formula
  a <- wsum[i]
  b <- -vwsum[i]
  d <- l
  delta <- -4*b^3*d - 27*a^2*d^2
  delta0 <- b^2
  delta1 <- 2*b^3 + 27*a^2*d
  delta2 <- -27*a^2*delta
  Ci <- (delta1+0i + sqrt(delta2+0i))^(1/3) / 2^(1/3)
  u <- c(1, complex(1,-1,-sqrt(3))/2, complex(1,-1,sqrt(3))/2)
  roots <- -(b + u*Ci + delta0/(Ci*u))/3/a
  realroots <- as.numeric(Re(roots[which(v.allequal(Im(roots), 0))]))     # only use the real part of those roots that are real, i.e. have imaginary part 0
  posind <- realroots > 0
  if(any(posind)) taus[i] <- min(realroots[posind]) else taus[i] <- -1    # only use solutions that are positive
 }
 tauind <- v > taus & taus >= c(v[2:p], 0)                                # which solutions lie within the corresponding interval?
 bestsol <- treshlnonex(u.orig, 1e42, w.orig, nn)
 if(any(tauind)){
  for(i in which(tauind)){
   soli <- treshlnonex(u.orig, taus[i], w.orig, nn)
   if(f(soli) < f(bestsol)) bestsol <- soli
  }
 }
 bestsol
})


###################################################################################################


##### proximal operator for category-specific predictors with global effects
## l2norm and generic soft thresholding function for global effects
treshg <- cmpfun(function(u, lambda, w, nn){
 l <- nrow(u)
 u <- u[1,]
 if(nn) u <- pmax(u, 0)
 w <- sqrt(ncol(u)) * w
 st <- 1 - lambda * w / l2normraw(u)
 u <- max(st, 0) * u
 matrix(rep(u, l), nrow=l, byrow=T)
})

 
###################################################################################################

#### proximal operator for category-specific predictors with category-specific effects
## l2norm and generic soft thresholding function 
treshcs <- cmpfun(function(u, lambda, w, nn){
 if(nn) u <- pmax(u, 0)
 w <- sqrt(length(c(u))) * w                                                
 st <- 1 - lambda * w / l2normraw(u)                                            
 u <- max(st, 0) * u
 u
})

## soft thresholding function for ordinary lasso penalties.
treshlcs <- cmpfun(function(u, lambda, w, nn){
 if(nn) u <- pmax(u, 0)
 w <- sqrt(ncol(u)) * w
 st <- 1 - lambda * w / sqrt(diag(tcrossprod(u)))
 u <- c(pmax(st, 0)) * u
 u
})
 
## sparse thresholding for sparse group lasso. input u must be a Q*P matrix.
sptreshcs <- cmpfun(function(u, lambda, lambdasp, w, wsp, nn){
 if(lambdasp > 0) u <- treshlrefx(u, lambdasp, wsp, nn)
 if(lambda > 0) u <- treshrefx(u, lambda, w, nn)
 u
})

###################################################################################################
## sparse thresholding for sparse group lasso combined with GFL based on ADMM.
spGFLtreshnonex <- cmpfun(function(v, lambdal2, lambdal1, lambdaF, wl2, wl1, wF, R, nn, trace=F){

 ## algorithmic tuning parameter
 rho <- 0.1 #1/lambdaF #0.1
 rhoscale <- 1 + 1/rho
 
 ## initialize stuff
 ## adjust weights for problem size
 wl2a <- sqrt(length(c(v))) * wl2
 wl1a <- sqrt(ncol(v)) * wl1
 wFa <- sqrt(ncol(v)) * sqrt(nrow(R)) * wF

 ## a function for computing the objective function of spGFLtreshnonex
 objective <- function(x){0.5*sum((x - v)^2) + lambdal2*wl2a*l2normraw(x) +
                          lambdal1*sum(wl1a*rowL2normraw(x)) + lambdaF*wFa*l2normraw(R%*%x)}

 ## algorithmic tuning parameters
 rel.tol <- 5e-3
 iter.max <- 200
 
 ## precompute and store t(R) and R'R
 tR <- t(R)
 #RtR <- crossprod(R)

 ## store the correct tresholding function belonging to the sparse group lasso part of the penalty
 if(lambdal2 * lambdal1 > 0) xtresh <- function(x) sptreshnonex(x, lambdal2/rhoscale, lambdal1/rhoscale, wl2, wl1, nn)
 else if(lambdal2 == 0) xtresh <- function(x) treshlnonex(x, lambdal1/rhoscale, wl1, nn)
 else if(lambdal1 == 0) xtresh <- function(x) treshnonex(x, lambdal2/rhoscale, wl2, nn)

 ## store the correct tresholding function for the z-update
 ztresh <- function(z) treshnonex(z, lambdaF, wF, nn)

 ## initialize the algorithmic quantities
 if(nn) v <- pmax(v, 0)
 x <- xtresh(v)
 z <- ztresh(R%*%x)
 u <- matrix(0, nrow=nrow(R), ncol=ncol(v))               # this is the scaled residual
 
 ## prepare for the loop
 iter.count <- 0
 d.fn <- 1e30
 d.fn.counter <- 0
 d.fn.converged <- F
 fn.val <- fn.val.old<- 1e30
 primal.converged <- dual.converged <- F

 #### the main ADMM loop ####
 while(!d.fn.converged || !primal.converged || !dual.converged){
  iter.count <- iter.count + 1
  if(iter.count > iter.max){
   if(trace){
    warning(paste("maximum number of iterations reached in spGFLtreshnonex"))
   }
   break
  }

  z.old <- z
  fn.val.old <- fn.val

  linpoint <- x  -  rho * tR %*% (R%*%x - z + u)
  xpoint <- (rho*v + linpoint)/(1 + rho)
  
  x <- xtresh(xpoint)                           # x-update
  Rx <- R%*%x
  
  zpoint <- Rx + u
  z <- ztresh(zpoint)                           # z-update

  primal.res <- Rx - z                          # primal residual
  u <- u + primal.res                           # u-update
  
  dual.res <- tR%*%(z - z.old)                  # dual residual                 #/ lambdaF ### hier eventuell durch lambdaF oder rho teilen?
  
  if(l2normraw(primal.res) <= rel.tol * max(c(l2normraw(Rx), l2normraw(z)))){
   primal.converged <- T
  }
  #if(l2normraw(dual.res) <= rel.tol * l2normraw(tR%*%u)){  ### hier eventuell noch rechte Seite mal lambdaF oder rho nehmen?
   dual.converged <- T
  #}

  if(!d.fn.converged){
   fn.val <- objective(x)                          # value of the objective
   d.fn <- abs(fn.val - fn.val.old)/(rel.tol/10 + abs(fn.val))
   if(d.fn < rel.tol) d.fn.converged <- T
  }
 } ## end while(!d.fn.converged ...)
 x
})

###################################################################################################
## sparse thresholding for sparse group lasso combined with L1 fusion based on ADMM.
spFGLtreshnonex <- cmpfun(function(v, lambdal2, lambdal1, lambdaF, wl2, wl1, wF, R, nn, trace=F){

 ## algorithmic tuning parameter
 rho <- 0.1 #1/lambdaF #0.1
 rhoscale <- 1 + 1/rho

 ## initialize stuff
 ## adjust weights for problem size
 wl2a <- sqrt(length(c(v))) * wl2
 wl1a <- sqrt(ncol(v)) * wl1
 wFa <- sqrt(ncol(v)) * wF

 ## a function for computing the objective function of spFGLtreshnonex
 objective <- function(x){0.5*sum((x - v)^2) + lambdal2*wl2a*l2normraw(x) +
                          lambdal1*sum(wl1a*rowL2normraw(x)) + lambdaF*sum(wFa*rowL2normraw(R%*%x))}

 ## algorithmic tuning parameters
 rel.tol <- 5e-3
 iter.max <- 200

 ## precompute and store t(R) and R'R
 tR <- t(R)
 #RtR <- crossprod(R)

 ## store the correct tresholding function belonging to the sparse group lasso part of the penalty
 if(lambdal2 * lambdal1 > 0) xtresh <- function(x) sptreshnonex(x, lambdal2/rhoscale, lambdal1/rhoscale, wl2, wl1, nn)
 else if(lambdal2 == 0) xtresh <- function(x) treshlnonex(x, lambdal1/rhoscale, wl1, nn)
 else if(lambdal1 == 0) xtresh <- function(x) treshnonex(x, lambdal2/rhoscale, wl2, nn)

 ## store the correct tresholding function for the z-update
 ztresh <- function(z) treshlnonex(z, lambdaF, wF, nn)

 ## initialize the algorithmic quantities
 if(nn) v <- pmax(v, 0)
 x <- xtresh(v)
 z <- ztresh(R%*%x)
 u <- matrix(0, nrow=nrow(R), ncol=ncol(v))               # this is the scaled residual

 ## prepare for the loop
 iter.count <- 0
 d.fn <- 1e30
 d.fn.counter <- 0
 d.fn.converged <- F
 fn.val <- fn.val.old<- 1e30
 primal.converged <- dual.converged <- F

 #### the main ADMM loop ####
 while(!d.fn.converged || !primal.converged || !dual.converged){
  iter.count <- iter.count + 1
  if(iter.count > iter.max){
   if(trace){
    warning(paste("maximum number of iterations reached in spGFLtreshnonex"))
   }
   break
  }

  z.old <- z
  fn.val.old <- fn.val

  linpoint <- x  -  rho * tR %*% (R%*%x - z + u)
  xpoint <- (rho*v + linpoint)/(1 + rho)

  x <- xtresh(xpoint)                           # x-update
  Rx <- R%*%x

  zpoint <- Rx + u
  z <- ztresh(zpoint)                           # z-update

  primal.res <- Rx - z                          # primal residual
  u <- u + primal.res                           # u-update

  dual.res <- tR%*%(z - z.old)                  # dual residual                 #/ lambdaF ### hier eventuell durch lambdaF oder rho teilen?

  if(l2normraw(primal.res) <= rel.tol * max(c(l2normraw(Rx), l2normraw(z)))){
   primal.converged <- T
  }
  #if(l2normraw(dual.res) <= rel.tol * l2normraw(tR%*%u)){  ### hier eventuell noch rechte Seite mal lambdaF oder rho nehmen?
   dual.converged <- T
  #}

  if(!d.fn.converged){
   fn.val <- objective(x)                          # value of the objective
   d.fn <- abs(fn.val - fn.val.old)/(rel.tol/10 + abs(fn.val))
   if(d.fn < rel.tol) d.fn.converged <- T
  }
 } ## end while(!d.fn.converged ...)
 x
})


##################################################################################################




###################################################################################################
## sparse thresholding for a "rowwise" L1 fusion based on ADMM.
rowFLtreshnonex <- cmpfun(function(v, lambdaF, wF, D, nn, trace=F){

 ## algorithmic tuning parameter
 rho <- 0.1 #1/lambdaF #0.1
 rhoscale <- 1 + 1/rho

 ## initialize stuff
 v <- t(v)              ## this enables the whole MRSP machinery, which is based on "columnwise" operations, for the rowwise case we consider here.
 tD <- t(D)
 
 ## a function for computing the objective function of rowFLtreshnonex
 objective <- function(x){0.5*sum((x - v)^2) + lambdaF*sum(wF*abs(D%*%x))}

 ## algorithmic tuning parameters
 rel.tol <- 5e-3
 iter.max <- 200

 ## store the correct tresholding function belonging to the sparse group lasso part of the penalty
 ##### here, this is not used. to make as few changes to existing code as possible, the following is used:
 xtresh <- function(x){ x }

 ## store the correct tresholding function for the z-update
 ztresh <- function(z) treshlnonex(z, lambdaF, wF, nn)

 ## initialize the algorithmic quantities
 if(nn) v <- pmax(v, 0)
 x <- xtresh(v)
 z <- ztresh(D%*%x)
 u <- matrix(0, nrow=nrow(D), ncol=ncol(v))               # this is the scaled residual

 ## prepare for the loop
 iter.count <- 0
 d.fn <- 1e30
 d.fn.counter <- 0
 d.fn.converged <- F
 fn.val <- fn.val.old<- 1e30
 primal.converged <- dual.converged <- F

 #### the main ADMM loop ####
 while(!d.fn.converged || !primal.converged || !dual.converged){
  iter.count <- iter.count + 1
  if(iter.count > iter.max){
   if(trace){
    warning(paste("maximum number of iterations reached in spGFLtreshnonex"))
   }
   break
  }

  z.old <- z
  fn.val.old <- fn.val

  linpoint <- x  -  rho * tD %*% (D%*%x - z + u)
  xpoint <- (rho*v + linpoint)/(1 + rho)

  x <- xtresh(xpoint)                           # x-update
  Dx <- D%*%x

  zpoint <- Dx + u
  z <- ztresh(zpoint)                           # z-update

  primal.res <- Dx - z                          # primal residual
  u <- u + primal.res                           # u-update

  dual.res <- tD%*%(z - z.old)                  # dual residual                 #/ lambdaF ### hier eventuell durch lambdaF oder rho teilen?

  if(l2normraw(primal.res) <= rel.tol * max(c(l2normraw(Dx), l2normraw(z)))){
   primal.converged <- T
  }
  #if(l2normraw(dual.res) <= rel.tol * l2normraw(tR%*%u)){  ### hier eventuell noch rechte Seite mal lambdaF oder rho nehmen?
   dual.converged <- T
  #}

  if(!d.fn.converged){
   fn.val <- objective(x)                          # value of the objective
   d.fn <- abs(fn.val - fn.val.old)/(rel.tol/10 + abs(fn.val))
   if(d.fn < rel.tol) d.fn.converged <- T
  }
 } ## end while(!d.fn.converged ...)
 t(x)  ## we have to transpose to revert the transposition from the beginning...
})


##################################################################################################



## helperfunction for compiling
fistaProximal.MRSP <- function(coef, tuning, penweights, penindex, grpindex,
                              Proximal.args=NULL, Proximal.control=NULL, ...)
{
 if(Proximal.args$constraint == "symmetric"){                                   
  tresh <- treshsymx
  treshl <- treshlsymx
  sptresh <- sptreshsymx
  Psplinetresh <- Psplinetreshsymx
  colPsplinetresh <- colPsplinetreshsymx
  spFGLtresh <- spFGLtreshsymx
 }else if(is.numeric(Proximal.args$constraint)){
  tresh <- treshrefx
  treshl <- treshlrefx
  sptresh <- sptreshrefx
  Psplinetresh <- Psplinetreshrefx
  colPsplinetresh <- colPsplinetreshrefx
 }else if(Proximal.args$constraint == "none"){
  tresh <- treshnonex
  treshl <- treshlnonex
  sptresh <- sptreshnonex
  spGFLtresh <- spGFLtreshnonex
  colPsplinetresh <- colPsplinetreshnonex
  spFGLtresh <- spFGLtreshnonex
  rowFLtresh <- rowFLtreshnonex
 } 
   
 if(Proximal.args$any.grouppen.x){
  for(j in Proximal.args$groups.grouppen.x){
   if(penweights[[1]][[1]][j] != 0){
    j.which <- which(grpindex[[1]] == j)
    coef[[1]][,j.which] <- tresh(coef[[1]][,j.which], tuning[[1]], penweights[[1]][[1]][j], Proximal.args$nonneg)
   }
  }
 }
 
 if(Proximal.args$any.spgrppen.x){
  for(j in Proximal.args$groups.spgrppen.x){
   if(penweights[[1]][[1]][j] != 0 | any(penweights[[2]][[1]][,j] != 0)){
    j.which <- which(grpindex[[1]] == j)
    coef[[1]][,j.which] <- sptresh(coef[[1]][,j.which, drop=F], tuning[[1]], tuning[[2]], penweights[[1]][[1]][j], penweights[[2]][[1]][,j], Proximal.args$nonneg)
   }
  }
 }
 
 if(Proximal.args$any.lassopen.x){
  for(j in Proximal.args$groups.lassopen.x){
   if(any(penweights[[2]][[1]][,j] != 0)){
    j.which <- which(grpindex[[1]] == j)
    coef[[1]][,j.which] <- treshl(coef[[1]][,j.which, drop=F], tuning[[2]], penweights[[2]][[1]][,j], Proximal.args$nonneg)
   }
  }
 }
 
 if(Proximal.args$any.spGFLpen.x){
  for(j in Proximal.args$groups.spGFLpen.x){
   #if(penweights[[1]][[1]][j] != 0 | any(penweights[[2]][[1]][,j] != 0) | penweights[[3]][[1]][j] != 0){
   if(penweights[[1]][[1]][j] != 0 | any(penweights[[2]][[1]][,j] != 0)){
    j.which <- which(grpindex[[1]] == j)
    if(tuning[[7]] > 0 && penweights[[3]][[1]][j] != 0){
     coef[[1]][,j.which] <- spGFLtresh(coef[[1]][,j.which], tuning[[1]], tuning[[2]], tuning[[7]], penweights[[1]][[1]][j],
                                       penweights[[2]][[1]][,j], penweights[[3]][[1]][j], Proximal.args$R, Proximal.args$nonneg)
    #coef[[1]][,j.which] <- sptresh(coef[[1]][,j.which, drop=F], tuning[[1]], tuning[[2]], penweights[[1]][[1]][j], penweights[[2]][[1]][,j])
    }else{
     coef[[1]][,j.which] <- sptresh(coef[[1]][,j.which, drop=F], tuning[[1]], tuning[[2]], penweights[[1]][[1]][j], penweights[[2]][[1]][,j], Proximal.args$nonneg)
    }
   }
  }
 }
 
 if(Proximal.args$any.Psplinepen.x){
  for(j in Proximal.args$groups.Psplinepen.x){
   j.which <- which(grpindex[[1]] == j)
   coef[[1]][,j.which] <- Psplinetresh(coef[[1]][,j.which], tuning[[7]], Proximal.args$Omega[[j]], Proximal.args$tD[[j]])
  }
 }
 
 if(Proximal.args$any.colPsplinepen.x){
  for(j in Proximal.args$groups.colPsplinepen.x){
   j.which <- which(grpindex[[1]] == j)
   reflevel <- Proximal.args$constraint
   if(is.numeric(reflevel)){
    coef[[1]][-reflevel, j.which] <- colPsplinetresh(coef[[1]][-reflevel, j.which], tuning[[7]], Proximal.args$R)
   }else{
    coef[[1]][,j.which] <- colPsplinetresh(coef[[1]][,j.which], tuning[[7]], Proximal.args$R)
   }
  }
 }
 
 if(Proximal.args$any.spFGLpen.x){
  for(j in Proximal.args$groups.spFGLpen.x){
   if(penweights[[1]][[1]][j] != 0 | any(penweights[[2]][[1]][,j] != 0)){
    j.which <- which(grpindex[[1]] == j)
    if(tuning[[7]] > 0 && any(penweights[[4]][[1]][,j] != 0)){
     coef[[1]][,j.which] <- spFGLtresh(coef[[1]][,j.which], tuning[[1]], tuning[[2]], tuning[[7]], penweights[[1]][[1]][j],
                                       penweights[[2]][[1]][,j], penweights[[4]][[1]][,j], Proximal.args$R, Proximal.args$nonneg)
    }else{
     coef[[1]][,j.which] <- sptresh(coef[[1]][,j.which, drop=F], tuning[[1]], tuning[[2]], penweights[[1]][[1]][j], penweights[[2]][[1]][,j], Proximal.args$nonneg)
    }
   }
  }
 }

 if(Proximal.args$any.rowFLpen.x){
  for(j in Proximal.args$groups.rowFLpen.x){
   if(any(penweights[[5]][[1]][[j]] != 0)){
    j.which <- which(grpindex[[1]] == j)
    coef[[1]][,j.which] <- rowFLtresh(coef[[1]][,j.which, drop=F], tuning[[7]], penweights[[5]][[1]][[j]], t(Proximal.args$tD[[j]]), Proximal.args$nonneg)
   }
  }
 }
  
 if(Proximal.args$any.ridgepen.x){
  coef[[1]][,Proximal.args$which.ridgepen.x] <- 1/(1 + tuning[[4]]) * coef[[1]][,Proximal.args$which.ridgepen.x] 
 }
 
 if(Proximal.args$constraint == "none" & Proximal.args$modelname %in% c("Multinomial Logit Model")){
  if(Proximal.args$any.notpen.x){
   coef[[1]][,Proximal.args$which.notpen.x] <- t(t(coef[[1]][,Proximal.args$which.notpen.x]) - treshmean(coef[[1]][,Proximal.args$which.notpen.x]))
  }
 }
   
 if(Proximal.args$hasV){
  lp <- length(penweights[[1]][[1]])

  if(Proximal.args$any.globalpen.V){
   for(j in Proximal.args$groups.globalpen.V){
    if(penweights[[1]][[2]][j - lp] != 0){
     j.which <- which(grpindex[[2]] == j)
     coef[[2]][,j.which] <- treshg(coef[[2]][,j.which], tuning[[3]], penweights[[1]][[2]][j - lp], Proximal.args$nonneg)
    }
   } 
  }
  
  if(Proximal.args$any.globalridgepen.V){
   coef[[2]][,Proximal.args$which.globalridgepen.V] <- 1/(1 + tuning[[4]]) * coef[[2]][,Proximal.args$which.globalridgepen.V] 
  }
  
  if(Proximal.args$any.catpen.V){
   for(j in Proximal.args$groups.catpen.V){
    if(penweights[[1]][[2]][j - lp] != 0){
     j.which <- which(grpindex[[2]] == j)
     coef[[2]][,j.which] <- treshcs(coef[[2]][,j.which], tuning[[3]], penweights[[1]][[2]][j - lp], Proximal.args$nonneg)
    } 
   }
  }
  
  if(Proximal.args$any.catspgrppen.V){
   for(j in Proximal.args$groups.catspgrppen.V){
    if(penweights[[1]][[2]][j - lp] != 0 | any(penweights[[2]][[2]][,j - lp] != 0)){
     j.which <- which(grpindex[[2]] == j)
     coef[[2]][,j.which] <- sptreshcs(coef[[2]][,j.which], tuning[[3]], tuning[[6]], penweights[[1]][[2]][j - lp], penweights[[2]][[2]][,j - lp], Proximal.args$nonneg)
    }
   }
  }   
  
  if(Proximal.args$any.catlassopen.V){
   for(j in Proximal.args$groups.catlassopen.V){
    if(any(penweights[[2]][[2]][,j - lp] != 0)){
     j.which <- which(grpindex[[2]] == j)
     coef[[2]][,j.which] <- treshlcs(coef[[2]][,j.which], tuning[[6]], penweights[[2]][[2]][,j - lp], Proximal.args$nonneg)
    }
   }
  }   
  
  if(Proximal.args$any.catridgepen.V){
   coef[[2]][,Proximal.args$which.catridgepen.V] <- 1/(1 + tuning[[4]]) * coef[[2]][,Proximal.args$which.catridgepen.V]
  } 
 } #end 'hasV'

 if(Proximal.args$isordinal){
  if(Proximal.args$any.globalpen.o){
   for(j in Proximal.args$groups.globalpen.o){
    if(penweights[[1]][[1]][j] != 0){
     j.which <- which(grpindex[[1]] == j)
     coef[[1]][,j.which] <- treshg(coef[[1]][,j.which], tuning[[1]], penweights[[1]][[1]][j], Proximal.args$nonneg)
    }
   }
  }

  if(Proximal.args$any.globalridgepen.o){
   coef[[1]][,Proximal.args$which.globalridgepen.o] <- 1/(1 + tuning[[4]]) * coef[[1]][,Proximal.args$which.globalridgepen.o]
  }
 } #end 'isordinal'

 if(Proximal.args$ridgestabil){
  coef <- lapply(coef, function(u){u/(1 + tuning[[5]])})
 }

 if(!is.null(Proximal.args$extrastabil)){
  coef[[1]][,which(Proximal.args$extrastabil[[1]])] <- 1/(1 + tuning[[5]]) * coef[[1]][,which(Proximal.args$extrastabil[[1]])]
  if(Proximal.args$hasV){
   coef[[2]][,which(Proximal.args$extrastabil[[2]])] <- 1/(1 + tuning[[5]]) * coef[[2]][,which(Proximal.args$extrastabil[[2]])]
  }
 }
 
 coef  
}

# compile it
fistaProximal.MRSP.cmpd <- cmpfun(fistaProximal.MRSP)

# the true method
setMethod("fistaProximal", signature(coef = "MRSP.coef"),
          function(coef, tuning, penweights, penindex, grpindex,
                   Proximal.args=NULL, Proximal.control=NULL, ...)
{
 fistaProximal.MRSP.cmpd(coef = coef, tuning = tuning, penweights = penweights,
              penindex = penindex, grpindex = grpindex, Proximal.args = Proximal.args,
              Proximal.control = Proximal.control, ...)
})



## penalty
# helperfunction for compiling
penalty.MRSP <- function(coef, tuning, penweights, penindex, grpindex,
                         Proximal.args=NULL, ...)
{
 coefsum <- 0
 
 if(Proximal.args$constraint == "symmetric" | is.numeric(Proximal.args$constraint)){
  l2normx <- l2norm.x
 }else if(Proximal.args$constraint == "none"){
  l2normx <- l2norm
 }  
 
 if(Proximal.args$any.grouppen.x){
  for(j in Proximal.args$groups.grouppen.x){
   j.which <- which(grpindex[[1]] == j)
   coefsum <- coefsum + tuning[[1]] * penweights[[1]][[1]][j] * l2normx(coef[[1]][,j.which])
  }
 }
 
 if(Proximal.args$any.spgrppen.x){
  for(j in Proximal.args$groups.spgrppen.x){
   j.which <- which(grpindex[[1]] == j)
   coefsum <- coefsum + penweights[[1]][[1]][j] * tuning[[1]] * l2normx(coef[[1]][,j.which]) + tuning[[2]] * sum(penweights[[2]][[1]][,j] * rowL2norm(coef[[1]][,j.which, drop=F]))
  }
 }
 
 if(Proximal.args$any.lassopen.x){
  for(j in Proximal.args$groups.lassopen.x){
   j.which <- which(grpindex[[1]] == j)
   coefsum <- coefsum + tuning[[2]] * sum(penweights[[2]][[1]][,j] * rowL2norm(coef[[1]][,j.which, drop=F]))
  }
 }

 if(Proximal.args$any.spGFLpen.x){
  for(j in Proximal.args$groups.spGFLpen.x){
   j.which <- which(grpindex[[1]] == j)
   #pGFL <- 1# length(Proximal.args$groups.spGFLpen.x)
   #dualoptj <- Proximal.args$dualopt.pen.x[[j]]
   coefsum <- coefsum + tuning[[1]] * penweights[[1]][[1]][j] * l2normx(coef[[1]][,j.which]) +
                        tuning[[2]] * sum(penweights[[2]][[1]][,j] * rowL2norm(coef[[1]][,j.which, drop=F])) +
                        sqrt(ncol(coef[[1]][,j.which])) * sqrt(nrow(Proximal.args$R)) * tuning[[7]] * penweights[[3]][[1]][j] * l2normraw(Proximal.args$R %*% coef[[1]][,j.which])
                        ## here, we must add the value of the smooth approximation of the GFL penalty when evaluated with the dual-optimal parameter values
                        #sqrt(ncol(coef[[1]][,j.which])) * max(sum(dualoptj %*% coef[[1]][,j.which]) - Proximal.args$SPG.eps/pGFL/2 * sum(dualoptj^2), 0)
  }
 }
 
 if(Proximal.args$any.Psplinepen.x){
  for(j in Proximal.args$groups.Psplinepen.x){
   j.which <- which(grpindex[[1]] == j)
   coefsum <- coefsum + tuning[[7]] * sum((coef[[1]][,j.which]%*%Proximal.args$tD[[j]])^2)
  }
 }
 
 if(Proximal.args$any.colPsplinepen.x){
  for(j in Proximal.args$groups.colPsplinepen.x){
   j.which <- which(grpindex[[1]] == j)
   reflevel <- Proximal.args$constraint
   if(is.numeric(reflevel)){
    coefsum <- coefsum + tuning[[7]] * sum((Proximal.args$R%*%coef[[1]][-reflevel, j.which])^2)
   }else{
    coefsum <- coefsum + tuning[[7]] * sum((Proximal.args$R%*%coef[[1]][,j.which])^2)
   }
  }
 }

 if(Proximal.args$any.spFGLpen.x){
  for(j in Proximal.args$groups.spFGLpen.x){
   j.which <- which(grpindex[[1]] == j)
   coefsum <- coefsum + tuning[[1]] * penweights[[1]][[1]][j] * l2normx(coef[[1]][,j.which]) +
                        tuning[[2]] * sum(penweights[[2]][[1]][,j] * rowL2norm(coef[[1]][,j.which, drop=F])) +
                        sqrt(ncol(coef[[1]][,j.which])) * tuning[[7]] * sum(penweights[[4]][[1]][,j] * rowL2normraw(Proximal.args$R %*% coef[[1]][,j.which]))
  }
 }

 if(Proximal.args$any.rowFLpen.x){
  for(j in Proximal.args$groups.rowFLpen.x){
   j.which <- which(grpindex[[1]] == j)
   coefsum <- coefsum + tuning[[7]] * sum(penweights[[5]][[1]][[j]] * abs(t(Proximal.args$tD[[j]])%*%t(coef[[1]][,j.which, drop=F])))
  }
 }

 if(Proximal.args$any.ridgepen.x){
  coefsum <- coefsum + tuning[[4]]/2 * sum(coef[[1]][,Proximal.args$which.ridgepen.x]^2)
 }
 
 if(Proximal.args$hasV){
  lp <- length(penweights[[1]][[1]])
 
  if(Proximal.args$any.globalpen.V){
   for(j in Proximal.args$groups.globalpen.V){
    j.which <- which(grpindex[[2]] == j)
    coefsum <- coefsum + penweights[[1]][[2]][j - lp] * tuning[[3]] * l2norm(coef[[2]][1,j.which])  ## only consider the first row!
   } 
  }
  
  if(Proximal.args$any.globalridgepen.V){
   coefsum <- coefsum + tuning[[4]]/2 * sum(coef[[2]][1,Proximal.args$which.globalridgepen.V]^2) ## same!
  }
  
  if(Proximal.args$any.catpen.V){
   for(j in Proximal.args$groups.catpen.V){
    j.which <- which(grpindex[[2]] == j)
    coefsum <- coefsum + tuning[[3]] * penweights[[1]][[2]][j - lp] * l2norm(coef[[2]][,j.which])
   }
  }  
  
  if(Proximal.args$any.catspgrppen.V){
   for(j in Proximal.args$groups.catspgrppen.V){
    j.which <- which(grpindex[[2]] == j)
    coefsum <- coefsum + tuning[[3]] * penweights[[1]][[2]][j - lp] * l2norm(coef[[2]][,j.which]) + tuning[[6]] * sum(penweights[[2]][[2]][,j - lp] * rowL2norm(coef[[2]][,j.which]))
   }
  }  
  
  if(Proximal.args$any.catlassopen.V){
   for(j in Proximal.args$groups.catlassopen.V){
    j.which <- which(grpindex[[2]] == j)
    coefsum <- coefsum + tuning[[6]] * sum(penweights[[2]][[2]][,j - lp] * rowL2norm(coef[[2]][,j.which]))
   }
  }
  
  if(Proximal.args$any.catridgepen.V){
   coefsum <- coefsum + tuning[[4]]/2 * sum(coef[[2]][,Proximal.args$which.catridgepen.V]^2)
  } 
 } # end 'hasV'

 if(Proximal.args$isordinal){
  if(Proximal.args$any.globalpen.o){
   for(j in Proximal.args$groups.globalpen.o){
    j.which <- which(grpindex[[1]] == j)
    coefsum <- coefsum + penweights[[1]][[1]][j] * tuning[[1]] * l2norm(coef[[1]][1,j.which])  ## only consider the first row!
   }
  }

  if(Proximal.args$any.globalridgepen.o){
   coefsum <- coefsum + tuning[[4]]/2 * sum(coef[[1]][1,Proximal.args$which.globalridgepen.o]^2) ## same!
  }
 } # end 'isordinal'
 
 if(Proximal.args$ridgestabil){
  coefsum <- coefsum + tuning[[5]] * sum((do.call("c", lapply(coef, "c")))^2)
 }

 if(!is.null(Proximal.args$extrastabil)){
  coefsum <- coefsum + tuning[[5]] * sum(c(coef[[1]][,which(Proximal.args$extrastabil[[1]])])^2)
  if(Proximal.args$hasV){
   coefsum <- coefsum + tuning[[5]] * sum(c(coef[[2]][,which(Proximal.args$extrastabil[[2]])])^2)
  }
 }

 as.numeric(coefsum)      
}

# compile it
penalty.MRSP.cmpd <- cmpfun(penalty.MRSP)

# the true method
setMethod("penalty", signature(coef = "MRSP.coef"),
          function(coef, tuning, penweights, penindex, grpindex,
                   Proximal.args=NULL, ...)
{
 penalty.MRSP.cmpd(coef = coef, tuning = tuning, penweights = penweights,
                   penindex = penindex, grpindex = grpindex,
                   Proximal.args = Proximal.args)
})


## lossapprox
setMethod("lossapprox", signature(coefdiff = "MRSP.coef"),
          function(lS, grad, coefdiff, Proximal.args, ...)
{
 indg <- Proximal.args$indg
 if(length(indg) > 0 & length(grad) > 1 & length(coefdiff) > 1){
  grad[[2]][-1,indg] <- 0
  coefdiff[[2]][-1,indg] <- 0
 }
 sl.indg <- Proximal.args$sl.indg
 if(length(sl.indg) > 0){
  grad[[1]][-1,sl.indg] <- 0
  coefdiff[[1]][-1,sl.indg] <- 0
 }
 lS + Reduce('+', Map(function(u,v){sum(u*v)}, grad, coefdiff))
})  

## proximterm
setMethod("proximterm", signature(coefdiff = "MRSP.coef"),
          function(coefdiff, Proximal.args, ...)
{
 indg <- Proximal.args$indg
 if(length(indg) > 0 & length(coefdiff) > 1){
  coefdiff[[2]][-1,indg] <- 0
 }
 sl.indg <- Proximal.args$sl.indg
 if(length(sl.indg) > 0){
  coefdiff[[1]][-1,sl.indg] <- 0
 }
 Reduce('+', lapply(coefdiff, function(u) sum(u^2)))
})


## update.eta
setMethod("updateEta", signature(coef = "MRSP.coef"),
          function(dat, coef, offset, weights, ...)
{
 eta <- tcrossprod(dat$x, coef[[1]])
 eta <- eta + offset
 if(!is.null(dat$V))
  #eta <- eta + matrix(do.call("rbind", dat$V) %*% coef[[2]], nrow = nrow(eta), byrow = T)
  eta <- eta + mapply("%*%", dat$V, as.list(as.data.frame(t(coef[[2]]))))    
 eta
})

## update.penweights      # currently not used
setMethod("updatePenweights", signature(coef = "MRSP.coef"),
          function(penweights, coef, weights, penindex, grpindex,
                   Proximal.args=NULL, ...)
{
penweights
})


## SPG function for MRSP
SPG.MRSP <- function(coef, tuning, penweights, grpindex, Proximal.args, doGrad=T, ...)
{
 ## initialize the output object and its 'grad' entry
 SPGlist <- list()
 if(doGrad){
  SPGlist$grad <- lapply(coef, function(u){if(is.matrix(u)){u <- matrix(0, nrow=nrow(u), ncol=ncol(u))}
                                           else if(is.vector(u)){u <- rep(0, length(u))}
                                           else stop("object coef had an unusable form in function 'SPG'")})
 }else{SPGlist$grad <- list()} 

 ## initialize the objects that store the dual-optimal variables for later use
 SPGlist$dualopt <- list()
 dualopt.pen.x <- list(); length(dualopt.pen.x) <- length(unique(grpindex[[1]]))

 ## for parameter groups with a penalty that uses SPG,
 ## compute the dual optimal solution and add the correct
 ## quantity to the respective parts of grad:
 if(Proximal.args$any.spGFLpen.x){
  R <- Proximal.args$R
  eps <- Proximal.args$SPG.eps
  lambdaF <- tuning[[7]]
  for(j in Proximal.args$groups.spGFLpen.x){
   j.which <- which(grpindex[[1]] == j)
   w <- penweights[[3]][[1]][j]
   pGFL <- 1 #length(Proximal.args$groups.spGFLpen.x)
   dualopt <- l2project(R%*%(coef[[1]][,j.which]), w, lambdaF, eps, pGFL)
   if(doGrad) SPGlist$grad[[1]][,j.which] <- lambdaF*w*(t(R)%*%dualopt)
   dualopt.pen.x[[j]] <- lambdaF*w*(t(dualopt)%*%R)
  }
 }
 SPGlist$dualopt$dualopt.pen.x <- dualopt.pen.x
 return(SPGlist)
}

## compile it
SPG.MRSP.cmpd <- cmpfun(SPG.MRSP)

## now set the actual SPG method
setMethod("SPG", signature(coef = "MRSP.coef"),
          function(coef, tuning, penweights, grpindex, Proximal.args, doGrad=T, ...)
{
 SPG.MRSP.cmpd(coef = coef, tuning = tuning, penweights = penweights,
               grpindex = grpindex, Proximal.args = Proximal.args, doGrad=doGrad)
})

## check if Smoothing Proximal Gradient has to be used
setMethod("useSPG", signature(coef = "MRSP.coef"),
          function(coef, penweights, Proximal.args, ...)
{
 F #if(any(c(Proximal.args$any.spGFLpen.x))) T else F
})

## compute the smoothed penalty if SPG is used
SPGsmoothpen.MRSP <- function(coef, dualopt, penweights, grpindex, tuning, Proximal.args, ...)
{
 pensum <- 0
 if(Proximal.args$any.spGFLpen.x){
  for(j in Proximal.args$groups.spGFLpen.x){
   j.which <- which(grpindex[[1]] == j)
   pGFL <- 1 # length(Proximal.args$groups.spGFLpen.x)
   dualoptj <- dualopt$dualopt.pen.x[[j]]
   pensum <- pensum + sqrt(ncol(coef[[1]][,j.which])) * tuning[[7]] * penweights[[3]][[1]][j] * l2normraw(Proximal.args$R %*% coef[[1]][,j.which])
   #pensum <- pensum + sqrt(ncol(coef[[1]][,j.which])) * max(sum(dualoptj %*% coef[[1]][,j.which]) - Proximal.args$SPG.eps/pGFL/2 * sum(dualoptj^2), 0)
  }
 }
 pensum
}

SPGsmoothpen.MRSP.cmpd <- cmpfun(SPGsmoothpen.MRSP)

setMethod("SPGsmoothpen", signature(coef = "MRSP.coef"),
          function(coef, dualopt, penweights, grpindex, tuning, Proximal.args, ...)
{
 SPGsmoothpen.MRSP.cmpd(coef=coef, dualopt=dualopt, penweights=penweights, grpindex=grpindex,
                        tuning=tuning, Proximal.args=Proximal.args)
})
