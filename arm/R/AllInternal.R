# some useful little functions
  #.round <- base:::round  
  sd.scalar <- function (x, ...) {sqrt(var(as.vector(x), ...))}
  wmean <- function (x, w, ...) {mean(x*w, ...)/mean(w, ...)}
  logit <- function (x) {log(x/(1-x))}
  .untriangle <- function (x) {x + t(x) - x*diag(nrow(as.matrix(x)))}




# new functions! 

as.matrix.VarCorr <- function (varc, useScale, digits){
# VarCorr function for lmer objects, altered as follows:
#   1.  specify rounding
#   2.  print statement at end is removed
#   3.  reMat is returned
#   4.  last line kept in reMat even when there's no error term
                  sc <- attr(varc, "sc")[[1]]   
                  if(is.na(sc)) sc <- 1 
#                  recorr <- lapply(varc, function(el) el@factors$correlation)
                  recorr <- lapply(varc, function(el) attr(el, "correlation"))
                  #reStdDev <- c(lapply(recorr, slot, "sd"), list(Residual = sc))
                  reStdDev <- c(lapply(varc, function(el) attr(el, "stddev")), list(Residual = sc))
                  reLens <- unlist(c(lapply(reStdDev, length)))
                  reMat <- array('', c(sum(reLens), 4),
                                 list(rep('', sum(reLens)),
                                      c("Groups", "Name", "Variance", "Std.Dev.")))
                  reMat[1+cumsum(reLens)-reLens, 1] <- names(reLens)
                  reMat[,2] <- c(unlist(lapply(reStdDev, names)), "")
#                  reMat[,3] <- format(unlist(reStdDev)^2, digits = digits)
#                  reMat[,4] <- format(unlist(reStdDev), digits = digits)
                  reMat[,3] <- fround(unlist(reStdDev)^2, digits)
                  reMat[,4] <- fround(unlist(reStdDev), digits)
                  if (any(reLens > 1)) {
                      maxlen <- max(reLens)
                      corr <-
                          do.call("rbind",
                                  lapply(recorr,
                                         function(x, maxlen) {
                                             x <- as(x, "matrix")
#                                             cc <- format(round(x, 3), nsmall = 3)
                                             cc <- fround (x, digits)
                                             cc[!lower.tri(cc)] <- ""
                                             nr <- dim(cc)[1]
                                             if (nr >= maxlen) return(cc)
                                             cbind(cc, matrix("", nr, maxlen-nr))
                                         }, maxlen))
                      colnames(corr) <- c("Corr", rep("", maxlen - 1))
                      reMat <- cbind(reMat, rbind(corr, rep("", ncol(corr))))
                  }
#                  if (!useScale) reMat <- reMat[-nrow(reMat),]
          if (useScale<0) reMat[nrow(reMat),] <- c ("No residual sd", rep("",ncol(reMat)-1))
          return (reMat)
      }


# rwish and dwish functions stolen from Martin and Quinn's MCMCpack

rwish <- function (v, S){
  if (!is.matrix(S)) 
        S <- matrix(S)
    if (nrow(S) != ncol(S)) {
        stop(message = "S not square in rwish().\n")
    }
    if (v < nrow(S)) {
        stop(message = "v is less than the dimension of S in rwish().\n")
    }
    p <- nrow(S)
    CC <- chol(S)
    Z <- matrix(0, p, p)
    diag(Z) <- sqrt(rchisq(p, v:(v - p + 1)))
    if (p > 1) {
        pseq <- 1:(p - 1)
        Z[rep(p * pseq, pseq) + unlist(lapply(pseq, seq))] <- rnorm(p * 
            (p - 1)/2)
    }
    return(crossprod(Z %*% CC))
}

dwish <- function (W, v, S) {
    if (!is.matrix(S)) 
        S <- matrix(S)
    if (nrow(S) != ncol(S)) {
        stop(message = "W not square in dwish()\n\n")
    }
    if (!is.matrix(W)) 
        S <- matrix(W)
    if (nrow(W) != ncol(W)) {
        stop(message = "W not square in dwish()\n\n")
    }
    if (nrow(S) != ncol(W)) {
        stop(message = "W and X of different dimensionality in dwish()\n\n")
    }
    if (v < nrow(S)) {
        stop(message = "v is less than the dimension of S in  dwish()\n\n")
    }
    k <- nrow(S)
    gammapart <- 1
    for (i in 1:k) {
        gammapart <- gammapart * gamma((v + 1 - i)/2)
    }
    denom <- gammapart * 2^(v * k/2) * pi^(k * (k - 1)/4)
    detS <- det(S)
    detW <- det(W)
    hold <- solve(S) %*% W
    tracehold <- sum(hold[row(hold) == col(hold)])
    num <- detS^(-v/2) * detW^((v - k - 1)/2) * exp(-1/2 * tracehold)
    return(num/denom)
}

# no visible binding~~~~~~~~~~~~~~~
# functions used to pass the check for bayespolr

pgumbel <- function(q, loc = 0, scale = 1, lower.tail = TRUE)
{
    q <- (q - loc)/scale
    p <- exp(-exp(-q))
    if (!lower.tail) 1 - p else p
}

dgumbel <- function (x, loc = 0, scale = 1, log = FALSE)
{
    d <- log(1/scale) - x - exp(-x)
    if (!log) exp(d) else d
}

# defin n to pass the bayesglm.fit and bayesglm.h.fit check
n <- NULL

# for mcplot
.pvalue <- function ( v1, v2 ){
  mean( ( sign( v1 - v2 ) + 1 ) / 2 )
}

.is.significant <- function ( p, alpha = 0.05 ){
  significant <- 0 + ( p > ( 1 - alpha ) ) - ( p < alpha )
  return( significant )
}


.weights.default <- function (object, ...) 
{
    wts <- object$weights
    if (is.null(wts)) 
        wts
    else napredict(object$na.action, wts)
}

#.sweep.inv <- function(G){
#  # sweeps a symmetric matrix on all positions
#  # (so inverts the matrix)
#  for(i in 1:nrow(G)) {
#    G <- .sweep.oper(G, i)
#  }
#  G
#}
#
#.sweep.oper <- function(G = theta, k = 1.){
#  # k is the sweep position
#  p <- dim(G)[1.]
#  H <- G
#  #first do generic elements (those that don't involve k)
#  H[] <- 0.
#  tmp <- matrix(G[, k], p, 1.) %*% matrix(G[, k], 1., p)
#  #now replace the row and col with index=k 
#  H <- G - tmp/G[k, k]
#  H[, k] <- G[, k]/G[k, k]
#  #now replace the (k,k) diagonal element 
#  H[k,  ] <- G[, k]/G[k, k]
#  # and we're done
#  H[k, k] <- -1./G[k, k]
#  H
#}
#
#
#.wls.all2 <- function(X, w = wts, Y = y, treat = Trt)
#{
#  #
#  # This produces coefficient estimates and both standard and robust variances 
#  # estimates for regression with weights
#  # the standard variance corresponds to a situation where an observation represents
#  # the mean of w observations
#  # the robust variance corresponds to a situation where weights represent 
#  # probability or sampling weights
#  #
#  # first put together the necessary data inputs
#  #
#  nunits <-  sum(w > 0)
#  k <-  ncol(X)
#  ## now the weights, properly normed
#  wn <-  w * (nunits/sum(w))
#  W <-  diag(wn * (nunits/sum(wn)))
#  #
#  # x prime x inverse (including weights)
#  vhat <-   - .sweep.inv((t(X) %*% W %*% X))
#  #
#  # estimated regression coefficients and variance for just the treatment coefficient
#  b <-  vhat %*% t(X) %*% W %*% Y
#  MSE <-  c(t(Y) %*% W %*% Y - t(b) %*% t(X) %*% W %*% Y)/(nunits - k)
#  var.std <-  (vhat * MSE)[2, 2]
#  #
#  ######  now for the robust variance calculations
#  # now a matrix where each row represents the contribution to the score
#  # for each observation
#  U <-  c((Y - X %*% b) * wn) * X
#  # finite sample adjustment
#  qc <-  nunits/(nunits - 2)
#  # the sum of outer products of each of the above score contributions for
#  # each person is calculated here
#  prodU <-  array(0, c(k, k, nunits))
#  for(i in 1:nunits) {
#    prodU[,  , i] <-  outer(U[i,  ], U[i,  ])
#  }
#  # putting it all together...
#  Vrob <-  qc * vhat %*% apply(prodU, c(1, 2), sum) %*% vhat
#  # and we pull off the variance just for the treatment effect 
#  var.rob <-  Vrob[2, 2]
#  ###############
#  results <-  c(var.std, var.rob, b[2])
#  results
#}
