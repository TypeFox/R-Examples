MLRC2 <- function(y, x, n.out=100, expand.grad=0.1, use.gam=FALSE, check.data=TRUE, lean=FALSE, n.cut=5, verbose=TRUE, ...)
{
  if (check.data) {
    if (any(apply(y, 1, sum) < 1.0E-8))
      stop(paste("Species data have zero abundances for the following rows:", paste(which(apply(y, 1, sum) < 1.0E-8), collapse=",")))
    if (any(apply(y, 2, sum) < 1.0E-8))
      stop(paste("Species data have zero abundances for the following columns:", paste(which(apply(y, 2, sum) < 1.0E-8), collapse=",")))
    if(n.cut < 5 & any(apply(y>0, 2, sum) < 5))
      warning("Trying to fit responses to some taxa with less than 5 occurrences - results may be unreliable")
  }
  if (any(y>1) | any (y<0))
    stop("Species data must be proportions between 0 and 1")
  fit <- MLRC2.fit(y=y, x=x, n.out=n.out, expand.grid=expand.grid, use.gam=use.gam, lean=lean, n.cut=n.cut, verbose=verbose, ...)
  xHat <- predict.internal.MLRC2(object=fit, y=y, lean=lean, ...) 
  call.print <- match.call()
  call.fit <- as.call(list(quote(MLRC2.fit), y=quote(y), x=quote(x), lean=FALSE))
  result <- c(fit, list(fitted.values=xHat, call.fit=call.fit, call.print=call.print, x=x))
  result$cv.summary <- list(cv.method="none")
  if (!lean) 
    result$y <- y
  class(result) <- "MLRC2" 
  result
}


MLRC2.fit <- function(y, x, n.out=100, expand.grad=0.1, use.gam=FALSE, check.data=TRUE, lean=FALSE, n.cut=5, verbose=TRUE, ...) {
#  require(mgcv)
  glr <- function(y, e, pred) {
    gfit <- glm(y ~ V1 + V2 + I(V2^2) + I(V2^2) + V1*V2, data=e, family = quasibinomial(link=logit), maxit = 100)
    p2 <- predict.glm(gfit, pred, type="response")
  }
  
  if (ncol(x) != 2)
    stop("x must contain 2 columns")
  nms <- colnames(x)
  colnames(x) <- c("V1", "V2")
  r.V1 <- range(x[, 1])
  exp.V1 <- (r.V1[2] - r.V1[1]) * expand.grad
  xpred1=seq(r.V1[1]-exp.V1, r.V1[2]+exp.V1, length.out=n.out)
  r.V2 <- range(x[, 2])
  exp.V2 <- (r.V2[2] - r.V2[1]) * expand.grad
  xpred2=seq(r.V2[1]-exp.V2, r.V2[2]+exp.V2, length.out=n.out)
  pred <- expand.grid(V1=xpred1, V2=xpred2)
#  res <- matrix(0, ncol = ncol(y), nrow = nrow(pred))
  
#  for (spec in 1:ncol(y))  {
#    if (use.gam) {
#      gfit <- gam(y[, spec] ~ s(V1) + s(V2), data=x, family = quasibinomial(link=logit), maxit = 100)
#      p2 <- predict.gam(gfit, pred, type="response")
#    }
#    else {
#      gfit <- glm(y[, spec] ~ V1 + V2 + I(V2^2) + I(V2^2) + V1*V2, data=x, family = quasibinomial(link=logit), maxit = #100)
#      p2 <- predict.glm(gfit, pred, type="response")
#    }
#    res[, spec] <- p2
#  }

  res <- apply(y, 2, glr, e=x, pred=pred)
  colnames(res) <- colnames(y)
  colnames(pred) <- nms
  list(spp=res, env=pred)
}

predict.internal.MLRC2 <- function(object, y, lean=FALSE, ...)
{
  if (!lean) {
    if (ncol(object$spp) != ncol(y))
      stop("Number of columns different in y and model in predict.internal.MLRC2")
  }
  nnn <- nrow(y)
  nn <- nrow(object$spp)
  LL.res <- matrix(nrow = nn, ncol = nnn)
  p <- log(object$spp)
  ppp <- log(1-object$spp)
  LL.res <- as.matrix(p) %*% t(y) + as.matrix(ppp) %*% t(1.0-y)
  LL.res[is.na(LL.res)] <- -1.0E30
  res3 <- object$env[apply(LL.res, 2, order, decreasing = TRUE)[1, ], 1:2]
  rownames(res3) <- rownames(y)
  colnames(res3) <- colnames(object$env)
  as.matrix(res3)
}

predict.MLRC2 <- function(object, newdata=NULL, sse=FALSE, nboot=100, match.data=TRUE, verbose=TRUE, ...) {
  if (!is.null(newdata))
    if (any(newdata < 0) | any(newdata > 1))
      stop("newdata must be proportions between 0 and 1")
  .predict(object=object, newdata=newdata, sse=sse, nboot=nboot, match.data=match.data, verbose=verbose, ...)
}

performance.MLRC2 <- function(object, ...) {
#  .performance(object, ...)
  RMSE <- apply(residuals(object), 2, .rmse)
  R2 <- diag(apply(object$fitted.values, 2, .r2, x=object$x))
  result <- cbind(RMSE, R2)
  out <- list(object=result)
  if (object$cv.summary$cv.method != "none") {
    RMSE <- apply(object$x - object$predicted, 2, .rmse)
    R2 <- diag(apply(object$predicted, 2, .r2, x=object$x))
    result.cv <- cbind(RMSE, R2)
    out$crossval <- result.cv
  }
  out
}

print.MLRC2 <- function(x, ...) 
{
  cat("\n")
  cat("Method : Maximum Likelihood using Response Curves \n")
  cat("Call   : ")
  cat(paste(deparse(x$call.print), "\n\n"))
  cat(paste("No. samples        :", length(x$x), "\n"))
  cat(paste("No. species        :", ncol(x$spp), "\n"))
  .print.crossval(x)
  cat("\nPerformance:\n")
  .print.performance(x)
  cat("\n")
}

summary.MLRC2 <- function(object, full=FALSE, ...) 
{
  print(object, ...)
  if (object$cv.summary$cv.method == "none")
    fitted <- as.data.frame(object$fitted.values)     
  else
    fitted <- as.data.frame(object$fitted.values, object$predicted)     
  cat("\nFitted values\n")
}

fitted.MLRC2 <- function(object, ...) {
  object$fitted.values
}

residuals.MLRC2 <- function(object, cv=FALSE, ...) {
  if (cv == FALSE)
    return (object$x - object$fitted.values)
  else {
    if (object$cv.summary$cv.method == "none")
      stop("Object does not contain cross validation results")
    return (object$residuals.cv)
  }
}

coef.MLRC2 <- function(object, ...) {
  n <- ncol(object$spp)
  res <- matrix(rep(NA, n))
  rownames(res) <- colnames(object$spp)
  return(res)
}
