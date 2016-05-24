ml.truncreg <- function(param, X, y, gradient = FALSE, hessian = FALSE, fit = FALSE, point, direction){
  beta <- param[1:ncol(X)]
  sigma <- param[length(param)]
  bX <- as.numeric(crossprod(beta,t(X)))
  resid <- y-bX
  if (direction == "left"){
    trunc <- bX-point
    sgn <- 1
  }
  else{
    trunc <- point-bX
    sgn <- -1
  }
  mills <- dnorm(trunc/sigma)/pnorm(trunc/sigma)

  lnL <- sum(log(dnorm(resid/sigma))-log(sigma)
               -log(pnorm(trunc/sigma)))
  if (gradient){
    gbX <- resid/sigma^2 -sgn/sigma*mills
    gsigma <- resid^2/sigma^3-1/sigma+trunc/sigma^2*mills
    gradi <- cbind(gbX*X,as.numeric(gsigma))
    attr(lnL,"gradient") <- gradi
  }
  if (fit){
    attr(lnL,"fit") <- bX
  }
  if (hessian){
    bb <- mills*(trunc/sigma+mills)/sigma^2-1/sigma^2
    ss <- -3*resid^2/sigma^4+1/sigma^2+trunc^2/sigma^4*mills*(mills+trunc/sigma)-
      2*trunc*mills/sigma^3
    bs <- -2*resid/sigma^3+sgn*(mills/sigma^2-trunc/sigma^3*mills*(mills+trunc/sigma))
    bb <- crossprod(bb*X,X)
    bs <- apply(bs*X,2,sum)
    ss <- sum(ss)
    h <- rbind(cbind(bb,bs),c(bs,ss))
    attr(lnL,"hessian") <- h
  }
  lnL
}

truncreg <- function(formula, data, subset, weights, na.action, point = 0, direction = "left",
  model = TRUE, y = FALSE, x = FALSE, ...)
{
  formula.type <- TRUE
  if (class(formula[[3]]) == "name"){
    X <- try(eval(formula[[3]],sys.frame(which=-3)),silent = TRUE)
    if (class(X) == "matrix") formula.type <- FALSE else formula.type <- TRUE
  }
  
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  
  if (formula.type){
    m <- match(c("formula", "data", "subset", "na.action","weights"),
               names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    X <- model.matrix(formula, data = mf)
    Y <- model.response(mf)
    mt <- attr(mf, "terms")
  }
  else{
    Y <- eval(formula[[2]], sys.frame(which=-3))
    mt <- terms(formula)
  }
  
  ## process options
  direction <- match.arg(direction, c("left", "right"))
  point <- rep(point, length.out = length(Y))
  
  result <- truncreg.fit(X, Y, point, direction, ...)
  result$call <- cl
  result$terms <- mt
  if(model) result$model <- mf
  if(y) result$y <- Y
  if(x) result$x <- X
  result
}

truncreg.fit <- function(X, y, point, direction, ...)
{
  ## check input
  if(direction == "left" & any(y < point)) stop("response not truncated, contains observations < 'point'")
  if(direction == "right" & any(y > point)) stop("response not truncated, contains observations > 'point'")

  dots <- list(...)
  if (is.null(dots$method)) method <- "bfgs" else method <- dots$method
  if (is.null(dots$iterlim)) iterlim <- 50 else iterlim <- dots$iterlim
  if (is.null(dots$print.level)) print.level <- 0 else print.level <- dots$print.level
  
  oldoptions <- options(warn=-1)
  on.exit(options(oldoptions))
  start.time <- proc.time()

  f <- function(param)  ml.truncreg(param,  X = X, y = y,
                                    gradient = FALSE, hessian = FALSE,
                                    fit = FALSE, point = point,
                                    direction = direction)
  g <- function(param){
    attr(ml.truncreg(param, X = X, y = y,
                     gradient = TRUE, hessian = FALSE,
                     fit = FALSE, point = point,
                     direction = direction),"gradient")
  } 
  h <- function(param){
    attr(ml.truncreg(param, X = X, y = y,
                     gradient = FALSE, hessian = TRUE,
                     fit = FALSE, point = point,
                     direction = direction),"hessian")
  } 

  linmod <- lm.fit(X, y)
  start <- c(linmod$coefficients, sqrt(sum(linmod$residuals^2)/linmod$df.residual))
  maxl <- maxLik(f, g, h, start = start, method = method,
                 iterlim = iterlim, print.level = print.level)
  grad.conv <- g(maxl$estimate)
  coefficients <- maxl$estimate
  vcov <- -solve(maxl$hessian)
  fit <- attr(ml.truncreg(coefficients, X = X, y = y,
                          gradient = FALSE, hessian = FALSE,
                          fit = TRUE, point = point,
                          direction = direction), "fit")
  names(fit) <- rownames(X)
  logLik <- maxl$maximum
  attr(logLik,"df") <- length(coefficients)
  hessian <- maxl$hessian
  convergence.OK <- maxl$code<=2
  elaps.time <- proc.time() - start.time
  nb.iter <- maxl$iterations
  eps <- maxl$gradient%*%solve(-maxl$hessian)%*%maxl$gradient

  est.stat <- list(elaps.time = elaps.time,
                   nb.iter = nb.iter,
                   eps = eps,
                   method = maxl$type,
                   message = maxl$message
                   )
  class(est.stat) <- "est.stat"
  coef.names <- c(colnames(X),"sigma")
  names(coefficients) <- rownames(vcov) <- colnames(vcov) <- coef.names

  result <- list(coefficients = coefficients,
                 vcov = vcov,
                 fitted.values = fit,
                 logLik = logLik,
                 gradient = grad.conv,
		 nobs = length(y),
                 call = NULL,
		 terms = NULL,
                 model = NULL,
		 y = NULL,
		 x = NULL,
		 point = if(isTRUE(all.equal(rep(point[1], length(point)), point))) point[1] else point,
		 direction = direction,
                 est.stat = est.stat
                 )
  class(result) <- c("truncreg", "maxLik")
  result
}

fitted.truncreg <- function(object, ...){
  object$fitted.values
}

residuals.truncreg <- function(object, ...){
  model.frame(object)[[1]]-fitted(object)
}

coef.truncreg <- function(object, ...){
  object$coefficients
}

vcov.truncreg <- function(object, ...){
  object$vcov
}

logLik.truncreg <- function(object, ...)
  structure(object$logLik, df = length(object$coefficients), nobs = object$nobs, class = "logLik")

print.truncreg <- function (x, digits = max(3, getOption("digits") - 2), width = getOption("width"), ...){
  cat("\nCall:\n", deparse(x$call), "\n\n", sep = "")
  if (length(coef(x))) {
    cat("Coefficients:\n")
    print.default(format(coef(x), digits = digits), print.gap = 2, 
                  quote = FALSE)
  }
  else cat("No coefficients\n")
  cat("\n")
  invisible(x)
}

summary.truncreg <- function (object,...){
  b <- coef(object)
  std.err <- sqrt(diag(vcov(object)))
  z <- b/std.err
  p <- 2*(1-pnorm(abs(z)))
  CoefTable <- cbind(b,std.err,z,p)
  colnames(CoefTable) <- c("Estimate","Std. Error","t-value","Pr(>|t|)")
  object$coefficients <- CoefTable
  class(object) <- c("summary.truncreg","truncreg")
  return(object)
}

print.summary.truncreg <- function(x,digits= max(3, getOption("digits") - 2),width=getOption("width"),...){
  cat("\nCall:\n")
  print(x$call)
  cat("\n")
  cat("\nCoefficients :\n")
  printCoefmat(x$coefficients,digits=digits)
  cat("\n")
  df <- attr(x$logLik,"df")
  cat(paste("Log-Likelihood: ",
            signif(x$logLik,digits),
            " on ",df," Df\n",sep=""))
  invisible(x)
}

model.frame.truncreg <- function(formula, ...) {
  if(!is.null(formula$model)) return(formula$model)
  NextMethod()
}

model.matrix.truncreg <- function(object, ...)
  if(!is.null(object$x)) object$x else model.matrix(object$terms, model.frame(object), ...)


predict.truncreg <- function(object, newdata = NULL, na.action = na.pass, ...) 
{
  if(missing(newdata)) {
    rval <- object$fitted.values
  } else {
    mt <- delete.response(object$terms)
    X <- model.matrix(mt, model.frame(mt, newdata, na.action = na.action))
    rval <- drop(X %*% head(object$coefficients, -1))
  }
  return(rval)
}
