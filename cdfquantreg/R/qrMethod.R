#' @title S3 Methods for getting output from fitted cdfqr Objects.
#' @description Give the Gradient Function for CDF-Quantile Distribution Modles
#' @param object The fitted cdfqr model.
#' @param x The fitted cdfqr model.
#' @param type The parts of coefficients or variance-covariance matrix to be extracted.
#' @param formula. Changes to the formula. See \code{\link[Formula]{update.Formula}} for details.
#' @param evaluate If true evaluate the new updated model else return the call for the new model.
#' @param digits Number of digits to be retained in printed output.
#' @param ... Pass onto other functions or currently ignored
#' @examples
#' data(cdfqrExampleData)
#' fit <- cdfquantreg(crc99 ~ vert | confl, 't2','t2', data = JurorData)
#' 
#' summary(fit)
#' print(fit)
#' logLik(fit)
#' coef(fit)
#' deviance(fit)
#' vcov(fit)
#' 
#' #Update the model
#' fit2 <- update(fit, crc99 ~ vert*confl | confl)
#' summary(fit2)
#' 
#' @method summary cdfqr
#' @export
summary.cdfqr <- function(object, ...) {
  print(object)
  #return(object)
}

#' @method print cdfqr
#' @export
#' @rdname summary.cdfqr
print.cdfqr <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  
  cat("Family: ", paste(x$family$fd, x$family$sd, collapse = "-"), "\n")
  cat("Call: ", deparse(x$call), "\n", fill = TRUE)
  
  if (length(x$coefficients$location)) {
    
    mutitile <- "Mu coefficients (Location submodel)\n"
    if (x$family$fd == 'km'){
      mutitile <- "1st shape parameter coefficients (log link)\n"
    }
    cat(mutitile)
    printCoefmat(x$coefficients$location, digits=digits,signif.legend = FALSE)
    cat("\n")
  } else{
    mutitile2 <- "No coefficients in location submodel\n\n"
    if (x$family$fd == 'km'){
      mutitile2 <- "No coefficient for the 1st shape parameter (log link)\n"
    }
    cat(mutitile2)
  }
  
  if (length(x$coefficients$dispersion)) {
    if (length(x$coefficients$dispersion)) {
      Sigmatitile <- "Sigma coefficients (Dispersion submodel)\n"
      if (x$family$fd == 'km'){
        Sigmatitile <- "2nd shape parameter coefficients (log link)\n"
      }
      cat(Sigmatitile)
      printCoefmat(x$coefficients$dispersion, digits=digits)
      cat("\n")
    } else{
      Sigmatitile2 <- "No coefficients in Dispersion submodel\n"
      if (x$family$fd == 'km'){
        Sigmatitile2 <- "No coefficient for the 2nd shape parameter (log link)\n"
      }
      cat(Sigmatitile2)
    }
  }
  
  cat("Converge: ", x$converged, fill = TRUE)
  cat("Log-Likelihood: ", round(x$logLik, digits=digits), "\n", fill = TRUE)
  cat("Gradient: ", round(x$grad, digits=digits), "\n", fill = TRUE)
  invisible(x)
  
}

#' @method logLik cdfqr
#' @export
#' @rdname summary.cdfqr
logLik.cdfqr <- function(object, ...) {
  object$logLik
}

#' @method deviance cdfqr
#' @export
#' @rdname summary.cdfqr
deviance.cdfqr <- function(object, ...) {
  object$deviance
}

#' @method coef cdfqr
#' @export
#' @rdname summary.cdfqr
coef.cdfqr <- function(object, type = c("full","mean","sigma"), ...) {
  type <- match.arg(type)
  mean <- object$coefficients$location[,1]
  sigma <- object$coefficients$dispersion[,1]
  
  coef <- switch(type, 
               full = {c(mean,sigma)}, 
               mean = {mean}, 
               sigma = {sigma})
  return(coef)
}

#' @method vcov cdfqr
#' @export
#' @rdname summary.cdfqr
vcov.cdfqr <- function(object, type = c("full","mean","sigma"), ...) {
  
  vc <- object$vcov
  k <- nrow(object$coefficients$location)
  m <- nrow(object$coefficients$dispersion)
  
  type = match.arg(type)
  
  mean <- vc[seq.int(length.out = k), seq.int(length.out = k), drop = FALSE]
  
  dispersion <- vc[seq.int(length.out = m) + k, seq.int(length.out = m) + k, drop = FALSE]
 
  vc <- switch(type, full = {object$vcov}, 
               mean = {mean}, 
              sigma = {dispersion})
  return(vc)
}

#' @method update cdfqr
#' @export
#' @rdname summary.cdfqr
update.cdfqr <- function(object, formula., ..., evaluate = TRUE) {
  call <- object$call
  if (!missing(formula.)) 
    call$formula <- formula(update(Formula::Formula(formula(object)), formula.))
  
  extras <- match.call(expand.dots = FALSE)$...

  if (length(extras)) {
    existing <- !is.na(match(names(extras), names(call)))
    for (a in names(extras)[existing]) call[[a]] <- extras[[a]]
    if (any(!existing)) {
      call <- c(as.list(call), extras[!existing])
      call <- as.call(call)
    }
  }
  if (evaluate) 
    eval(call, parent.frame())
  else call
} 

#' @method formula cdfqr
#' @export
#' @rdname summary.cdfqr
formula.cdfqr <- function(x, ...) {
  call <- x$call
  return(formula(call$formula))
} 

