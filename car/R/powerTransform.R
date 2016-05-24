# revision history
# 2009-09-16: added ... argument to print.summary.powerTransform. J. Fox
# 2015-02-02: added 'gamma' argument to get transformation of (U + gamma)
# 2015-08-10: added estimateTransform as a generic function
# 2015-08-24: made 'family' an explicit argument to powerTransformation to clairfy man page.

### Power families:
basicPower <- function(U,lambda, gamma=NULL) {
 if(!is.null(gamma)) basicPower(t(t(as.matrix(U) + gamma)), lambda) else{
 bp1 <- function(U,lambda){
  if(any(U[!is.na(U)] <= 0)) stop("First argument must be strictly positive.")
  if (abs(lambda) <= 1.e-6) log(U) else (U^lambda)
  }
  out <- U
  out <- if(is.matrix(out) | is.data.frame(out)){
    if(is.null(colnames(out))) colnames(out) <- 
        paste("Z", 1:dim(out)[2],sep="")
    for (j in 1:ncol(out)) {out[, j] <- bp1(out[, j],lambda[j]) 
       colnames(out)[j] <- if(abs(lambda[j]) <= 1.e-6)
           paste("log(", colnames(out)[j],")", sep="") else
           paste(colnames(out)[j], round(lambda[j], 2), sep="^")}
    out}  else
    bp1(out, lambda)
  out}}
  
bcPower <- function(U, lambda, jacobian.adjusted=FALSE, gamma=NULL) {
 if(!is.null(gamma)) bcPower(t(t(as.matrix(U) + gamma)), lambda, jacobian.adjusted) else{
 bc1 <- function(U, lambda){
  if(any(U[!is.na(U)] <= 0)) stop("First argument must be strictly positive.")
  z <- if (abs(lambda) <= 1.e-6) log(U) else ((U^lambda) - 1)/lambda
  if (jacobian.adjusted == TRUE) {
    z * (exp(mean(log(U), na.rm=TRUE)))^(1-lambda)} else z
  }
  out <- U
  out <- if(is.matrix(out) | is.data.frame(out)){
    if(is.null(colnames(out))) colnames(out) <- 
        paste("Z", 1:dim(out)[2], sep="")
    for (j in 1:ncol(out)) {out[, j] <- bc1(out[, j], lambda[j]) }
    colnames(out) <- paste(colnames(out), round(lambda, 2), sep="^")
    out}  else
    bc1(out, lambda)
  out}}
  
yjPower <- function(U, lambda, jacobian.adjusted=FALSE) {
 yj1 <- function(U, lambda){
  nonnegs <- U >= 0
  z <- rep(NA, length(U))
  z[which(nonnegs)] <- bcPower(U[which(nonnegs)]+1, lambda, jacobian.adjusted=FALSE)
  z[which(!nonnegs)] <- -bcPower(-U[which(!nonnegs)]+1, 2-lambda, jacobian.adjusted=FALSE)
  if (jacobian.adjusted == TRUE)
        z * (exp(mean(log((1 + abs(U))^(2 * nonnegs - 1)), na.rm=TRUE)))^(1 -
            lambda)
    else z
  }
  out <- U
  out <- if(is.matrix(out) | is.data.frame(out)){
    if(is.null(colnames(out))) colnames(out) <- 
        paste("Z", 1:dim(out)[2], sep="")
    for (j in 1:ncol(out)) {out[, j] <- yj1(out[, j], lambda[j]) }
    colnames(out) <- paste(colnames(out), round(lambda, 2), sep="^")
    out}  else
    yj1(out, lambda)
  out}
  
powerTransform <- function(object, ...) UseMethod("powerTransform")

powerTransform.default <- function(object, family="bcPower", ...) {
   y <- object
   if(!inherits(y, "matrix") & !inherits(y, "data.frame")) {
       y <- matrix(y,ncol=1)
       colnames(y) <- c(paste(deparse(substitute(object))))}
   y <- na.omit(y)
   x <- rep(1, dim(y)[1]) 
   estimateTransform(x, y, NULL, family=family, ...)
   }                                    

powerTransform.lm <- function(object, family="bcPower", ...) {
    mf <- if(is.null(object$model)) 
            update(object, model=TRUE, method="model.frame")$model 
            else object$model
    mt <- attr(mf, "terms")
        y <- model.response(mf, "numeric")
    w <- as.vector(model.weights(mf))
    if (is.null(w)) w <- rep(1, dim(mf)[1])
    if (is.empty.model(mt)) {
        x <- matrix(rep(1,dim(mf)[1]), ncol=1) }
    else {
        x <- model.matrix(mt, mf, contrasts)   } 
  estimateTransform(x, y, w, family=family, ...)
  }                                                 
  
powerTransform.formula <- function(object, data, subset, weights, na.action, family="bcPower",
  ...) {
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("object", "data", "subset", "weights", "na.action"), 
         names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- as.name("model.frame")
    names(mf)[which(names(mf)=="object")] <- "formula"
    mf <- eval(mf, parent.frame())
    mt <- attr(mf, "terms")
    y <- model.response(mf, "numeric")
    w <- as.vector(model.weights(mf))
    if (is.null(w)) w <- rep(1, dim(mf)[1])
    if (is.empty.model(mt)) {
        x <- matrix(rep(1, dim(mf)[1]), ncol=1) }
    else {
        x <- model.matrix(mt, mf)   } 
  estimateTransform(x, y, w, family=family, ...)
  } 

estimateTransform <- function(X, Y, weights=NULL,
                              family="bcPower", start=NULL, method="L-BFGS-B", ...) {
  if(family == "skewPower") 
    estimateTransform.skewPower(X, Y, weights,  ...) else
      estimateTransform.default(X, Y, weights, family, start, method, ...)
}

# estimateTransform.default is renamed 'estimateTransform
estimateTransform.default <- function(X, Y, weights=NULL,
                                      family="bcPower", start=NULL, method="L-BFGS-B", ...) {
  fam <- match.fun(family)
  Y <- as.matrix(Y) # coerces Y to be a matrix.
  X <- as.matrix(X) # coerces X to be a matrix. 
  w <- if(is.null(weights)) 1 else sqrt(weights)
  nc <- dim(Y)[2]
  nr <- nrow(Y)
  xqr <- qr(w * X)       
  llik <- function(lambda){
    (nr/2)*log(((nr - 1)/nr) *
                 det(var(qr.resid(xqr, w*fam(Y, lambda, j=TRUE, ...)))))
  }
  llik1d <- function(lambda,Y){
    (nr/2)*log(((nr - 1)/nr) * var(qr.resid(xqr, w*fam(Y, lambda, j=TRUE, ...))))
  }
  if (is.null(start)) {
    start <- rep(1, nc)
    for (j in 1:nc){
      res<- suppressWarnings(optimize(
        f = function(lambda) llik1d(lambda,Y[ , j, drop=FALSE]),
        lower=-3, upper=+3))
      start[j] <- res$minimum
    }
  }
  res<-optim(start, llik, hessian=TRUE, method=method,  ...)
  if(res$convergence != 0)
    warning(paste("Convergence failure: return code =", res$convergence))
  res$start<-start
  res$lambda <- res$par
  names(res$lambda) <- 
    if (is.null(colnames(Y))) paste("Y", 1:dim(Y)[2], sep="")
  else colnames(Y)
  roundlam <- res$lambda
  stderr <- sqrt(diag(solve(res$hessian)))
  lamL <- roundlam - 1.96 * stderr
  lamU <- roundlam + 1.96 * stderr
  for (val in rev(c(1, 0, -1, .5, .33, -.5, -.33, 2, -2))) { 
    sel <- lamL <= val & val <= lamU 
    roundlam[sel] <- val
  }
  res$roundlam <- roundlam
  res$par <- NULL
  res$family<-family
  res$xqr <- xqr
  res$y <- Y
  res$x <- as.matrix(X)    
  res$weights <- weights
  class(res) <- "powerTransform"
  res
}
   
testTransform <- function(object, lambda) UseMethod("testTransform")   
   
testTransform.powerTransform <- function(object, lambda=rep(1, dim(object$y)[2])){
   fam <- match.fun(object$family)
   Y <- cbind(object$y) # coerces Y to be a matrix.
   nc <- dim(Y)[2]
   nr <- nrow(Y)
   lam <- if(length(lambda)==1) rep(lambda, nc) else lambda
   xqr <- object$xqr 
   w <- if(is.null(object$weights)) 1 else sqrt(object$weights)  
   llik <- function(lambda){
        (nr/2) * log(((nr - 1)/nr) * 
             det(var(qr.resid(xqr, w * fam(Y, lam, jacobian.adjusted=TRUE)))))
        }
   LR <- 2 * (llik(lambda) - object$value)
   df <- length(object$lambda)
   pval <- 1-pchisq(LR, df)
   out <- data.frame(LRT=LR, df=df, pval=pval)
   rownames(out) <- 
     c(paste("LR test, lambda = (",
             paste(round(lam, 2), collapse=" "), ")", sep=""))
   out}   

print.powerTransform<-function(x, ...) {
   cat("Estimated transformation parameters \n")
   print(x$lambda)
   invisible(x)}
      
summary.powerTransform<-function(object,...){
    one <- 1==length(object$lambda)
    label <- paste(object$family, 
       (if(one) "Transformation to Normality" else 
                "Transformations to Multinormality"), "\n\n")
    lambda<-object$lambda
    stderr<-sqrt(diag(solve(object$hessian)))
    df<-length(lambda) 
    result <- cbind(lambda, stderr, lambda - 1.96*stderr, lambda + 1.96*stderr)
    rownames(result)<-names(object$lambda)
    colnames(result)<-c("Est.Power", "Std.Err.", "Wald Lower Bound",
                        "Wald Upper Bound")
    tests <- testTransform(object, 0)
    tests <- rbind(tests, testTransform(object, 1))
    if ( !(all(object$roundlam==0) | all(object$roundlam==1) | 
        length(object$roundlam)==1 ))
           tests <- rbind(tests, testTransform(object, object$roundlam))
    out <-  list(label=label, result=result, tests=tests)
    class(out) <- "summary.powerTransform"
    out
    }
    
print.summary.powerTransform <- function(x,digits=4, ...) {
    cat(x$label)
    print(round(x$result, digits))
    cat("\nLikelihood ratio tests about transformation parameters\n")
    print(x$tests)
    }

coef.powerTransform <- function(object, round=FALSE, ...) 
  if(round==TRUE) object$roundlam else object$lambda
  
vcov.powerTransform <- function(object,...) {
  ans <- solve(object$hessian)
  rownames(ans) <- names(coef(object))
  colnames(ans) <- names(coef(object))
  ans}
  
plot.powerTransform <- function(x, z=NULL, round=TRUE, plot=pairs, ...){
  y <- match.fun(x$family)(x$y, coef(x, round=round))
  if (is.null(z)) plot(y, ...) else
   if (is.matrix(z) | is.data.frame(z)) plot(cbind(y, z), ...) else {
    y <- cbind(y,z)
    colnames(y)[dim(y)[2]] <- deparse(substitute(z))
    plot(y, ...) }
    }




  
  

                              