print.aftgee <- function(x, ...) {
  cat("Call:\n")
  print(x$call)
  cat("Coefficients:\n")
  print(x$coef.res)
  cat("\n Initial Estimator:\n")
  print(x$coef.init)
}

print.aftsrr <- function(x, ...) {
  cat("Call:\n")
  print(x$call)
  cat("Coefficients:\n")
  print(x$beta)
}

summary.aftgee <- function(object,...){
  z <- object
  if (class(z) != "aftgee"){
    stop("Most be aftgee class")
  }
  ans <- z["call"]
  TAB.ini <- NULL
  ## aftgee part
  est.gee <- z$coef.res
  se.gee <- sqrt(diag(z$var.res))
  est.temp.gee <- ifelse(se.gee == "NaN", "NaN", est.gee)
  z.val.gee <- as.numeric(est.temp.gee)/as.numeric(se.gee)
  TAB <- cbind(Estimate = round(est.gee, 3),
               StdErr = round(se.gee, 3),
               z.value = round(z.val.gee, 3),
               p.value = round(2 * pnorm(-abs(z.val.gee)), 3))
  rownames(TAB) <- names(z$coef.res)
  ## binit part
  est.ini <- z$coef.init
  res <- list(call=object$call, coefficients=TAB, binit = z$binit, iniEst = z$iniEst, est.ini = z$coef.init)
  class(res) <- "summary.aftgee"
  res
}

summary.aftsrr <- function(object,...){
  z <- object
  if (class(z) != "aftsrr"){
    stop("Most be aftsrr class")
  }
  ans <- z["call"]
  var.meth <- z$var.meth[z$var.meth %in% c("MB", "ZLCF", "ZLMB", "sHCF", "sHMB", "ISCF", "ISMB", "js")]
  se.count <- length(var.meth)
  se.name <- match(var.meth, names(z$covmat))
  ## if (z$intercept == TRUE) {
  ##     z$beta <- z$beta[-1]
  ##     for (i in se.name) {
  ##         z$covmat[[i]] <- z$covmat[[i]][-1, -1]
  ##     }
  ##     z$vari.name <- z$vari.name[-1]
  ## }
  ## se.covmat <- z$covmat[[se.name]]
  est.srr <- z$beta
  p <- length(z$beta)
  TAB.srr <- NULL
  ##  se.covmat <- list(NULL)
  ##  se.covmat[se.count + 1] <- NULL
  for (i in 1:se.count) {
      se.srr <- NA
      if (z$B != 0) {
          se.srr <- sqrt(diag(z$covmat[[se.name[i]]]))
      }
      z.val.srr <- as.numeric(est.srr)/as.numeric(se.srr)
      temp.srr <- cbind(Estimate = round(est.srr, 3), StdErr = round(se.srr, 3), z.value = round(z.val.srr, 3), p.value = round(2 * pnorm(-abs(z.val.srr)), 3))
      rownames(temp.srr) <- z$vari.name
      TAB.srr <- append(TAB.srr, list(temp.srr))
  }
  res <- list(call = object$call, coefficients = TAB.srr, var.name = names(z$covmat)[se.name])
  class(res) <- "summary.aftsrr"
  res
}

print.summary.aftgee <- function(x, ...){
  cat("Call:\n")
  print(x$call)
  cat("\n")
  cat("AFTGEE Estimator")
  cat("\n")
  printCoefmat(as.matrix(x$coefficients), P.values = TRUE, has.Pvalue = TRUE)
  ## if (is.numeric(x$binit) != TRUE) {
  ##     if (x$binit == "lm") {
  ##         cat("Initial Estimator from lm:")
  ##         cat("\n")
  ##         cat(format(round(as.numeric(x$est.ini), digits = 5), nsmall = 5))
  ##     }
  ##     if (x$binit == "srrgehan") {
  ##         cat("Initial Estimator from aftsrr with Gehan's weight:")
  ##         cat("\n")
  ##         cat(format(round(as.numeric(x$est.ini), digits = 5), nsmall = 5))
  ##     }
  ##     cat("\n")
  ##     cat("AFTGEE Estimator:")
  ##     cat("\n")
  ##     printCoefmat(x$coefficients, P.values = TRUE, has.Pvalue = TRUE)
  ## }
  ## if (is.numeric(x$binit) == TRUE) {
  ##     cat("Gehan Estimator:")
  ##     cat("\n")
  ##     cat(format(round(as.numeric(x$est.ini), digits = 5), nsmall = 5))
  ##     cat("\n")
  ##     cat("AFTGEE Estimator")
  ##     cat("\n")
  ##     printCoefmat(as.matrix(x$coefficients), P.values = TRUE, has.Pvalue = TRUE)
  ## }
}


print.summary.aftsrr <- function(x, ...){
  se.count <- length(x$var.name)
  cat("Call:\n")
  print(x$call)
  for (i in 1:se.count){
      cat("\n")
      cat("Variance Estimator:", as.character(x$var.name[i]))
      cat("\n")
      printCoefmat(as.data.frame(x$coefficients[i]), P.values = TRUE, has.Pvalue = TRUE)
  }
}


coef.aftsrr <- function(object, ...){
  z <- object
  if (class(z) != "aftsrr"){
    stop("Most be aftsrr class")
  }
  ans <- z["call"]
  out <- z$beta
  names(out) <- z$vari.name
  out
}


residuals.aftsrr <- function(object, ...){
  z <- object
  if (class(z) != "aftsrr"){
    stop("Most be aftsrr class")
  }
  ans <- z["call"]
  out <- log(z$y[,1]) - z$x %*% z$beta
  out
}


vcov.aftsrr <- function(object, ...){
  z <- object
  if (class(z) != "aftsrr"){
    stop("Most be aftsrr class")
  }
  ans <- z["call"]
  var.meth <- z$var.meth[z$var.meth %in% c("MB", "ZLCF", "ZLMB", "sHCF", "sHMB", "ISCF", "ISMB", "js")]
  se.count <- length(var.meth)
  se.name <- match(var.meth, names(z$covmat))
  p <- length(z$beta)
  TAB.srr <- NULL
  out <- list(NULL)
  out[se.count + 1] <- NULL
  names(out) <- z$var.meth
  for (i in 1:se.count) {
      se.srr <- z$covmat[[se.name[i]]]
      rownames(se.srr) <- z$vari.name
      colnames(se.srr) <- z$vari.name
      out[[i]] <- se.srr
  }
  out
}

coef.aftgee <- function(object, ...){
  z <- object
  if (class(z) != "aftgee"){
    stop("Most be aftgee class")
  }
  ans <- z["call"]
  out <- z$coef.res
  out
}


vcov.aftgee <- function(object, ...){
  z <- object
  if (class(z) != "aftgee"){
    stop("Most be aftgee class")
  }
  ans <- z["call"]
  out <- z$var.res
  out
}


residuals.aftgee <- function(object, ...){
  z <- object
  if (class(z) != "aftgee"){
    stop("Most be aftgee class")
  }
  ans <- z["call"]
  out <- log(z$y) - z$x %*% z$coef.res
  out
}


predict.aftsrr <- function(object, newdata = NULL, se.fit = FALSE, type = "lp", ...){
  z <- object
  out <- NULL
  if (is.null(newdata)) {
      out$fit <- as.numeric(z$x %*% z$beta)
      if (type == "response") {
          out$fit <- as.numeric(exp(out$fit))
      }
  }

  if (!is.null(newdata)) {
      n <- as.matrix(newdata, ncol = length(z$beta))
      out$fit <- as.numeric(n %*% z$beta)
      if (type == "response") {
          out$fit <- as.numeric(exp(out$fit))
      }
  }
  if (se.fit == TRUE) {
      var.meth <- z$var.meth[z$var.meth %in% c("MB", "ZLCF", "ZLMB", "sHCF", "sHMB", "ISCF", "ISMB", "js")]
      se.count <- length(var.meth)
      se.name <- match(var.meth, names(z$covmat))
      p <- length(z$beta)
      TAB.srr <- NULL
      var <- list(NULL)
      var[se.count + 1] <- NULL
          names(var) <- z$var.meth
      for (i in 1:se.count) {
          se.srr <- z$covmat[[se.name[i]]]
          ## rownames(se.srr) <- z$vari.name
          ## colnames(se.srr) <- z$vari.name
          if (is.null(newdata)) {
              var[[i]] <- as.numeric(sqrt(diag(z$x %*% se.srr %*% t(z$x))))
          }
          if (!is.null(newdata)) {
              var[[i]] <- as.numeric(sqrt(diag(n %*% se.srr %*% t(n))))
          }
      }
      out$se.fit <- var
      if (type == "response") {
          out$se.fit <- lapply(out$se.fit, function(x) out$fit * x)
      }
  }
  out
}

predict.aftgee <- function(object, newdata = NULL, se.fit = FALSE, ...){
    z <- object
    out <- NULL
    if (class(z) != "aftgee"){
        stop("Most be aftgee class")
    }
    ans <- z["call"]
    if (is.null(newdata)) {
        out$fit <- z$x %*% z$coef.res
    }

    if (!is.null(newdata)) {
        n <- as.matrix(newdata, ncol = length(z$coef.res))
        if (z$intercept == TRUE & ncol(n) < length(z$coef.res)) {
            n <- cbind(1, n)
        }
        out$fit <- n %*% z$coef.res
        if (se.fit == TRUE) {
            out$se.fit <- sqrt(diag(n %*% z$var.res %*% t(n)))
        }
    }

    out
}
