# ------------------------------------
# Author: Andreas Alfons
#         Erasmus University Rotterdam
# ------------------------------------

#' @method print bicSelect
#' @export
print.bicSelect <- function(x, best = TRUE, ...) {
  # print BIC values
  cat("\nBIC:\n")
  print(x$values, ...)
  # print optimal model if requested
  if(isTRUE(best)) {
    best <- x$best
    text <- "Index of best model:"
    if(length(best) == 1) {
      best <- as.matrix(best)
      dimnames(best) <- list(text, "")
    } else cat("\n", text, "\n", sep="")
    print(best, ...)
  }
  # return object invisibly
  invisible(x)
}

#' @method print fitSelect
#' @export
print.fitSelect <- function(x, ...) {
  # indices of the best reweighted and raw fit
  cat("Index of best model:\n")
  print(x$best, ...)
  # return object invisibly
  invisible(x)
}

#' @method print perrySeqModel
#' @export
#' @import perry
print.perrySeqModel <- function(x, ...) {
  # print prediction error results
  perry:::print.perrySelect(x, best=FALSE, ...)
  # print optimal step
  sOpt <- as.matrix(fits(x)[x$best])
  cat(sprintf("\nOptimal step: %d\n", sOpt))
  # print final model
  cat("\nFinal model:\n")
  print(x$finalModel, ...)
  # return object invisibly
  invisible(x)
}

#' @method print perrySparseLTS
#' @export
#' @import perry
print.perrySparseLTS <- function(x, ...) {
  # print prediction error results
  perry:::print.perryTuning(x, best=FALSE, final=FALSE, ...)
  # print optimal value for penalty parameter
  optimalLambda <- x$tuning[x$best, "lambda"]
  names(optimalLambda) <- names(x$best)
  cat("\nOptimal lambda:\n")
  print(optimalLambda, ...)
  # print final model
  cat("\nFinal model:\n")
  print(x$finalModel, ...)
  # return object invisibly
  invisible(x)
}

#' @method print seqModel
#' @export
print.seqModel <- function(x, zeros = FALSE, best = TRUE, ...) {
  # print function call
  if(!is.null(call <- x$call)) {
    cat("\nCall:\n")
    dput(x$call)
  }
  # print predictor sequence
  active <- t(x$active)
  steps <- seq_len(ncol(active))
  text <- if(inherits(x, "grplars")) "Group" else "Var"
  dimnames(active) <- list(text, steps)
  cat("\nSequence of moves:\n")
  print(active, ...)
  # print coefficients of optimal submodel
  sOpt <- getSOpt(x)
  if(is.null(sOpt)) {
    sOpt <- x$s  # only one step
    text <- c("Coefficients:", "Step:")
  } else text <- c("Coefficients of optimal submodel:", "Optimal step:")
  cat("\n", text[1], "\n", sep="")
  print(coef(x, zeros=zeros), ...)
  # print optimal step
  if(isTRUE(best)) cat("\n", text[2], sprintf(" %d\n", sOpt), sep="")
  # return object invisibly
  invisible(x)
}

#' @method print sparseLTS
#' @export
print.sparseLTS <- function(x, fit = c("reweighted", "raw", "both"),
                            zeros = FALSE, ...) {
  # initializations
  fit <- match.arg(fit)
  lambda <- x$lambda
  sOpt <- getSOpt(x, fit=fit)
  # print function call
  if(!is.null(call <- x$call)) {
    cat("\nCall:\n")
    dput(x$call)
  }
  # print coefficients of sparse LTS model
  coefficients <- coef(x, fit=fit, zeros=zeros)
  if(length(lambda) == 1) cat("\nCoefficients:\n")
  else {
    if(fit == "both") colnames(coefficients) <- c("reweighted", "raw")
    cat("\nCoefficients of optimal submodel:\n")
  }
  print(coefficients, ...)
  # print penalty parameter and robust scale estimate
  scale <- getScale(x, fit=fit)
  if(length(lambda) == 1) {
    text <- c("Penalty parameter:", "Residual scale estimate:")
    if(fit == "both") {
      lambda <- as.matrix(lambda)
      dimnames(lambda) <- list(text[1], "")
      print(lambda, ...)
      scale <- t(scale)
      dimnames(scale) <- list(text[2], colnames(coefficients))
      print(scale, ...)
    } else {
      info <- matrix(c(lambda, scale), 2, 1)
      dimnames(info) <- list(text, "")
      print(info, ...)
    }
  } else {
    lambda <- lambda[sOpt]
    info <- rbind(lambda, scale)
    text <- c("Optimal penalty parameter:", "Residual scale estimate:")
    if(fit == "both") {
      dimnames(info) <- list(text, colnames(coefficients))
      cat("\n")
    } else dimnames(info) <- list(text, "")
    print(info, ...)
  }
  # return object invisibly
  invisible(x)
}

#' @method print tslarsP
#' @export
print.tslarsP <- function(x, ...) {
  # print "grplars" model
  print.seqModel(x, best=FALSE, ...)
  # print optimal step and lag length
  sOpt <- getSOpt(x)
  if(is.null(sOpt)) {
    sOpt <- x$s  # only one step
    text <- "Step:"
  } else text <- "Optimal step:"
  info <- rbind(sOpt, x$p)
  dimnames(info) <- list(c(text, "Lag length:"), "")
  print(info, ...)
  # return object invisibly
  invisible(x)
}

#' @method print tslars
#' @export
print.tslars <- function(x, ...) {
  # print "grplars" model with optimal lag length
  pOpt <- x$pOpt
  xOpt <- x$pFit[[pOpt]]
  xOpt$call <- x$call
  print.seqModel(xOpt, best = FALSE, ...)
  # print optimal step and lag length
  sOpt <- getSOpt(xOpt)
  if(is.null(sOpt)) {
    sOpt <- xOpt$s  # only one step
    text <- "Step:"
  } else text <- "Optimal step:"
  info <- rbind(sOpt, pOpt)
  dimnames(info) <- list(c(text, "Optimal lag length:"), "")
  print(info, ...)
  # return object invisibly
  invisible(x)
}
