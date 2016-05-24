# Feb 11 2015.  skewPower family of transformations
#  skewPower:  Evaluate transformation at (lambda, gamma);
#    Jacobian.adjustment optional
#    U Univeriate or multivariate.
#  skewPowerllik:  Evaluate the skew log-likelihood at (lambda, gamma)
# April 15, 2015, corrected handing of invHess submatrices in estimateTransform.skewPower    


## ------------------------------------------------------------------------
skewPower <- function(U, lambda, jacobian.adjusted=FALSE, gamma) { 
  if(is.matrix(U)){
    if(dim(U)[2] != length(lambda) | dim(U)[2] != length(gamma)) 
      stop("gamma and lambda must have length equal to number of columns in U")
  } else {
    if(length(gamma) != 1 | length(lambda) != 1) 
      stop("gamma and lambda must be length 1")
  } 
  if(any(gamma < 0)) stop("gamma must be >= 0")
  hc1 <- function(U, lambda, gamma){
    if(abs(gamma) <= 1.e-10 & any(U[!is.na(U)] <= 0)) 
      stop("First argument must be strictly positive if gamma = 0.")
    s <- sqrt(U^2 + gamma^2)
    z <- if (abs(lambda) <= 1.e-10) 
              log(.5*(U + s)) else ((.5*(U + s))^lambda - 1)/lambda
    if (jacobian.adjusted == TRUE) {
      Jn <- (.5^lambda) * 
        (exp((lambda - 1) * mean(log(U + s), na.rm=TRUE))) *
        (exp(mean(log(1 + U/s), na.rm=TRUE)))
      z <- z/Jn}
    z
  }
  out <- U
  out <- if(is.matrix(out) | is.data.frame(out)){
    if(is.null(colnames(out))) 
       colnames(out) <- paste("Z", 1:dim(out)[2], sep="")
    for (j in 1:ncol(out)) {out[, j] <- hc1(out[, j], lambda[j], gamma[j]) }
    colnames(out) <- paste(colnames(out), "(",round(lambda, 2), ",",round(gamma, 1),")", sep="")
#    colnames(out) <- paste(colnames(out), round(lambda, 2), sep="^")
    out}  else
      hc1(out, lambda, gamma)
  out}

## Evaluate skew llik at (lambda, gamma)-----------------------------------
skewPowerllik <- function(X, Y, weights=NULL, lambda, gamma, xqr=NULL) {
  Y <- as.matrix(Y) # coerces Y to be a matrix.
  w <- if(is.null(weights)) 1 else sqrt(weights)  
  xqr <- if(is.null(xqr)){qr(w * as.matrix(X))} else xqr
  nr <- nrow(Y)
  f <- -(nr/2)*log(((nr - 1)/nr) *
      det(as.matrix(var(qr.resid(xqr, w * skewPower(Y, lambda, jacobian.adjusted=TRUE, gamma=gamma))))))
  list(lambda=lambda, gamma=gamma, llik=f)
}

## maximize skewPowerllik for fixed gamma--------------------------------------
skewPowerllikprofile.lambda <- function(X, Y, weights=NULL, gamma, xqr=NULL){
  xqr <- if(is.null(xqr)){
    w <- if(is.null(weights)) 1 else sqrt(weights)
    qr(w * as.matrix(X))
  } else xqr
  fn <- function(lam) skewPowerllik(NULL, Y, weights, lambda=lam, gamma=gamma, xqr=xqr)$llik
  if(dim(as.matrix(Y))[2] ==1){  
    f <- optimize(f=fn, interval=c(-3, 3), maximum=TRUE)
    list(lambda=f$maximum, gamma=gamma, llik=f$objective, invHess=solve(-optimHess(f$maximum, fn)))
  } else {
# get starting values
    par <- rep(0, length(gamma))
    for (j in 1:dim(Y)[2]) par[j] <- 
      skewPowerllikprofile.lambda(NULL, Y[, j, drop=FALSE], weights, gamma=gamma[j], xqr)$lambda
# optimize
    f <- optim(par=par, fn=fn, method="L-BFGS-B",
                  lower=rep(-3, length(gamma)), upper=rep(3, length(gamma)),
                  control=list(fnscale=-1))
    list(lambda=f$par, gamma=gamma, llik=f$value, invHess = solve(-optimHess(f$par, fn))) 
  }  
}

##  maximize skewPowerllik for fixed lambda-------------------------------------
skewPowerllikprofile.gamma <- function(X, Y, weights=NULL, lambda, xqr=NULL){ 
  xqr <- if(is.null(xqr)){
    w <- if(is.null(weights)) 1 else sqrt(weights)
    qr(w * as.matrix(X))
  } else xqr
 fn <- function(gam) skewPowerllik(NULL, Y, weights, lambda=lambda, gamma=gam, xqr=xqr)$llik 
 if(dim(as.matrix(Y))[2] == 1L){ 
   f <- optimize(f=fn, interval=c(0.01, max(Y)), maximum=TRUE)
   list(lambda=lambda, gamma=f$maximum, llik=f$objective, 
        invHess=solve(-optimHess(f$maximum, fn))) 
 } else {
# get starting values
   par <- rep(0, length(lambda))
   for (j in 1:dim(Y)[2])par[j] <- 
      skewPowerllikprofile.gamma(NULL, Y[, j, drop=FALSE], weights, lambda=lambda[j], xqr)$gamma
# optimize
   f <- optim(par=par, fn=fn, method="L-BFGS-B",
              lower=rep(.001, length(gamma)), upper=rep(Inf, length(gamma)),
              control=list(fnscale=-1))
   list(lambda=lambda, gamma=f$par, llik=f$value,
        invHess = solve(-optimHess(f$par, fn))) 
 }  
}

## skew lme, multivariate -----------------------------------------------------------
skewmle <- function(X, Y, weights=NULL, lambda=c(-3, 3), control=list(maxit=1000)) {
# Get a starting value for lambda by removing all negative rows 
# and doing Box-Cox Power:
  Y <- as.matrix(Y)
  sel <- apply(Y, 1, function(x) all(x > 0))
  weights <- if(is.null(weights)) rep(1, length(sel)) else weights
  lambda.start <- estimateTransform(as.matrix(X[sel, ]), as.matrix(Y[sel,]), weights[sel])$lambda
# Get a starting value for gamma using profiling
  xqr <- qr(sqrt(weights) * as.matrix(X))
  nc <- length(lambda.start)
  gamma.start <- skewPowerllikprofile.gamma(NULL, Y, weights, lambda.start, xqr=xqr)$gamma
# fnscale must be negative to maximize
  control$fnscale <- if(is.null(control$fnscale)) -1 else -abs(control$fnscale)
  fn <- function(param){
    lam <- param[1:nc]
    gam <- param[(nc+1):(2*nc)]
    skewPowerllik(NULL, Y, weights, lam, gam, xqr=xqr)$llik
  } 
  res <- optim(c(lambda.start, gamma.start), fn, method="L-BFGS-B", 
               lower=c(rep(lambda[1], nc), rep(.01, nc)), 
               upper=c(rep(lambda[2], nc), rep(+Inf, nc)), 
               control=control)
  if(res$convergence != 0)
    warning(paste("Convergence failure:", res$message))
  res$invHess <- solve(-optimHess(res$par, fn))
  res$ylabs <-  
    if (is.null(colnames(Y))) paste("Y", 1:dim(Y)[2], sep="") else colnames(Y)
  names(res$par) <- rownames(res$invHess) <- 
    colnames(res$invHess) <- c(paste("lam",res$ylabs, sep="."), paste("gam", res$ylabs, sep="."))
  res$xqr <- xqr
  res$y <- Y
  res$x <- X
  res$weights <- weights
  res
}

# Generic function is in powerTransform.R
estimateTransform.skewPower <- function(X, Y, weights=NULL, ...){
  res <- skewmle(X, Y, weights, ...)
  nc <- dim(as.matrix(Y))[2]
  res$lambda <- res$par[1:nc]
  res$gamma <- res$par[(nc+1):(2*nc)]
  roundlam <- res$lambda
  stderr <- sqrt(diag(res$invHess[1:nc, 1:nc, drop=FALSE]))
  stderr.gam <- sqrt(diag(res$invHess[(nc+1):(2*nc), (nc+1):(2*nc), drop=FALSE]))
  lamL <- roundlam - 1.96 * stderr
  lamU <- roundlam + 1.96 * stderr
  for (val in rev(c(1, 0, -1, .5, .33, -.5, -.33, 2, -2))) { 
    sel <- lamL <= val & val <= lamU 
    roundlam[sel] <- val
  }
  res$roundlam <- roundlam
  res$family <- "skewpowerTransform"
  class(res) <- c("skewpowerTransform", "powerTransform")
  res  
}

# generic is in powerTransform.R
testTransform.skewpowerTransform <- function(object, 
                      lambda=rep(1, dim(object$y)[2])){
  nc <- length(object$lambda)
  lam <- if(length(lambda)==1) rep(lambda, nc) else lambda  
  val <- skewPowerllikprofile.gamma(NULL, object$y, object$weights, lam, xqr=object$xqr)
  LR <- max(0, -2 * (val$llik - object$value))
  df <- nc
  pval <- 1-pchisq(LR, df)
  out <- data.frame(LRT=LR, df=df, pval=pval)
  rownames(out) <- 
    c(paste("LR test, lambda = (",
            paste(round(lam, 2), collapse=" "), ")", sep=""))
  out}   

# methods for skewpowerTransform objects

print.skewpowerTransform<-function(x, ...) {
  cat("Estimated transformation power, lambda\n")
  print(x$lambda)
  cat("Estimated transformation location, gamma\n")
  print(x$gamma)
  invisible(x)}

summary.skewpowerTransform<-function(object,...){
  nc <- length(object$lambda)
  label <- paste(if(nc==1) "Skew Power transformation to Normality" else 
                   "Skew Power Transformations to Multinormality", "\n\n")
  lambda <- object$lambda
  gamma <- object$gamma
  stderr <- sqrt(diag(object$invHess))
  stderr.gamma <- stderr[(nc+1):(2*nc)]
  stderr <- stderr[1:nc]
  result <- cbind(lambda, stderr, lambda - 1.96*stderr, lambda + 1.96*stderr)
  result.gamma <- cbind(gamma, stderr.gamma, pmax(gamma - 1.96*stderr.gamma, 0), gamma + 1.96*stderr.gamma)
  rownames(result) <- rownames(result.gamma) <- object$ylabs
  colnames(result) <- colnames(result.gamma) <- 
     c("Est.Power", "Std.Err.", "Wald Lower Bound", "Wald Upper Bound")
  colnames(result.gamma) <- 
    c("Est.gamma", "Std.Err.", "Wald Lower Bound", "Wald Upper Bound")
  tests <- testTransform(object, 0)
  tests <- rbind(tests, testTransform(object, 1))  
  if ( !(all(object$roundlam==0) | all(object$roundlam==1) | 
           length(object$roundlam)==1 ))
    tests <- rbind(tests, testTransform(object, object$roundlam))
  out <-  list(label=label, result=result, result.gamma=result.gamma, tests=tests)
  class(out) <- "summary.skewpowerTransform"
  out
}

print.summary.skewpowerTransform <- function(x,digits=4, ...) {
  cat(x$label)
  cat("\nEstimated power, lambda\n")
  print(round(x$result, digits))
  cat("\nEstimated location, gamma\n")
  print(round(x$result.gamma, digits)) 
  cat("\nLikelihood ratio tests about transformation parameters\n")
  print(x$tests) 
}

coef.skewpowerTransform <- function(object, param=c("both", "lambda", "gamma"), round=FALSE, ...){
  param <- match.arg(param)
  co <- cbind(if(round==TRUE) object$roundlam else object$lambda, object$gamma)
  dimnames(co) <- list(object$ylabs, c("lambda", "gamma"))
  switch(param, lambda = co[, 1], gamma=co[, 2], both= co)
}


vcov.skewpowerTransform <- function(object, param=c("both", "lambda", "gamma"), ...) {
  param <- match.arg(param)  
  nc <- length(object$lambda)
  switch(param, lambda=object$invHess[1:nc, 1:nc], gamma=object$invHess[(nc+1):(2*nc), (nc+1):(2*nc)],
         both=object$invHess)
}

plot.skewpowerTransform <- function(x, z=NULL, round=TRUE, plot=pairs, ...){
  y <- skewPower(x$y, lambda=coef(x, param="lambda"), 
                 jacobian.adjusted=FALSE, gamma=coef(x, param="gamma"))
  if (is.null(z)) plot(y, ...) else
    if (is.matrix(z) | is.data.frame(z)) plot(cbind(y, z), ...) else {
      y <- cbind(y, z)
      colnames(y)[dim(y)[2]] <- deparse(substitute(z))
      plot(y, ...) }
}

## ------------------------------------------------------------------------
contour.skewpowerTransform <- function(x, ksds=4, levels=c(.5, .95, .99, .999), 
              main="Skew Power Log-likelihood", ...){
  object <- x
  if(dim(object$y)[2] != 1L) stop("This function is for univariate Y only")
  q <- object$value - qchisq(levels, 2)/2
  se <- sqrt(diag(object$invHess))
  center <- c(object$lambda, object$gamma)
  x1 <- seq(object$lambda - ksds*se[1], object$lambda + ksds*se[1], length=100) 
  y <- seq(max(.01, object$gamma - ksds*se[2]), object$gamma + ksds*se[2], length=100)
  z <- matrix(0, nrow=length(x1), ncol=length(y))
  for (i in 1:length(x1)){
    for (j in 1:length(y)){
      z[i,j] <- skewPowerllik(NULL, object$y, object$weights, x1[i], y[j], xqr=object$xqr)$llik
    }
  }
  contour(x1, y, z, xlab=expression(lambda), ylab=expression(gamma), main=main, 
          nlevels=length(levels), levels=q, ...)
  points(center[1], center[2], pch=16, cex=1.25)
  text(center[1], center[2], as.character(round(object$value, 2)), pos=4, cex=.75)
}






