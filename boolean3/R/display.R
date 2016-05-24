### ----------------------------------------------------------------------------
### This file is part of boolean3
###
### Copyright (C) 2011--2014 Jason W. Morgan <morgan.746@osu.edu>
###
### boolean3 represents a substantial re-write of the original boolean package
### developed by Bear Braumoeller, Ben Goodrich, and Jacob Kline. This version
### was developed under the direction of Bear Braumoeller and with support from
### The Ohio State University's College of Social and Behavioral Sciences.
###
### boolean3 and is free software: you can redistribute it and/or modify it under
### the terms of the GNU General Public License as published by the Free Software
### Foundation, either version 3 of the License, or (at your option) any later
### version.
###
### This program is distributed in the hope that it will be useful, but WITHOUT
### ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
### FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
### details.
###
### You should have received a copy of the GNU General Public License along with
### this program.  If not, see <http://www.gnu.org/licenses/>.
###
### ----------------------------------------------------------------------------


##' Print a brief summary of a boolean model fit.
##'
##' This function prints a brief summary of the model fit to the user. The
##' details of a model fit are provided by \code{summary.boolean}. If the model
##' was not fit, this simply reports some summary statistics about the unfit
##' model.
##'
##' @method print boolean
##' @param x boolean object.
##' @param ... Additional parameters to be passed.
##' @return print brief summary of the model and the model fit (if available).
##' @author Jason W. Morgan (\email{morgan.746@@osu.edu})
##' @export print.boolean
print.boolean <- function(x, ...) {

  if (x$method[1] == "mcmc") {
    fit <- x$model.fit$mcmc
    par <- data.frame("mcmc"=fit$coef)
    rownames(par) <- x$coef.labels
    results <- rbind(par, "Chains"=fit$chains, "Thin"=fit$thin,
                     "Iterations"=fit$iterations)
  }
  else {
    par <- sapply(x$model.fit, function(mod) mod$par)
    rownames(par) <- x$coef.labels
    errids <- sapply(x$model.fit, function(mod) mod$conv)
    errmsg <- sapply(x$model.fit, function(mod) mod$message)
    lik <- do.call(c, unclass(logLik(x)))
    aic <- AIC(x)
    bic <- BIC(x)
    results <- rbind(par, "Log-likelihood"=lik, "AIC"=aic, "BIC"=bic)
  }

  cat("\n", paste(rep("=", 74), collapse=""), "\n", sep="")
  cat("Boolean model estimate\n\n")
  cat("Model specification:\n", as.character(x$call[2]), "\n")

  for (m in names(x$model$submods))
    cat("", m, "~", paste(x$model$submods[[m]], collapse=" + "), "\n")

  cat("\n")
  cat(paste(format(c("Observations", "Parameters")),
            paste(formatC(c(x$N, x$k), digits=3), "\n")), sep="")
  cat("\n")
  
  print(format(results, digits=3, width=4), quote=FALSE)

  cat("\n")
  if (!(x$method[1] == "mcmc")) {
    for (e in names(errids))
      if (errids[[e]] != 0) conv_msg(errids[[e]], msg=errmsg[[e]], method=e)
  }
  cat(paste(rep("=", 74), collapse=""), "\n\n") 
}


##' Create a summary object from a boolean estimate.
##'
##' This function calculates the standard model fit summary information given
##' the model object. 
##'
##' @method summary boolean
##' @param object object of class \code{boolean}.
##' @param ... Additional parameters to be passed.
##' @return A list of class \code{\link{boolsum-class}}.
##' @author Jason W. Morgan (\email{morgan.746@@osu.edu})
##' @export summary.boolean
summary.boolean <- function(object, ...) {

  ## Check that model.fit exists.
  if (is.null(object$model.fit))
      stop("'object' does not appear to contain a fit model", .call=FALSE)

  smry <- list(est.smry=c("Observations"=object$N, "Parameters"=object$k))
  method <- list()
  
  if (object$method[1] == "mcmc") {
    s <- summary(object$model.fit$mcmc$mcmc)
    ci <- apply(object$model.fit$mcmc$mcmc, 2, quantile, probs=c(0.025, 0.975))
    Coef <- cbind(s$statistics[,-3], t(ci))
    colnames(Coef) <- c("Mean", "SD", "SE(1)", "2.5%", "97.5%")
    rownames(Coef) <- object$coef.labels
    Smry <- c("Iterations"=object$model.fit$mcmc$iterations,
              "Thin"=object$model.fit$mcmc$thin,
              "Chains"=object$model.fit$mcmc$chains)
    method[["mcmc"]] <- list(coef=Coef, method.smry=Smry)
  }
  else {
    est <- coef(object)
    rownames(est) <- object$coef.labels
    se  <- sapply(vcov(object), function(x) sqrt(diag(x)))
    z   <- vapply(1:ncol(est), function(j) est[,j] / se[,j], rep(0, nrow(est)))
    pval <- vapply(1:ncol(z), function(j) pnorm(-abs(z[,j])), rep(0, nrow(z)))
    llik <- logLik(object)
    aic  <- AIC(object)
    bic  <- BIC(object)

    for (i in 1:ncol(est)) {
      Coef <- data.frame(Coef=est[,i], SE=se[,i], z=z[,i], "p-val"=pval[,i])
      Smry <- c("Log-likelihood"=as.numeric(llik[i]), "AIC"=aic[i],
                "BIC"=bic[i])
      method[[colnames(est)[i]]] <- list(coef=Coef, method.smry=Smry)
    }
  }
  smry[["methods"]] <- method
  
  ## Return only selected methods.
  methods <- list(...)$methods
  if (!is.null(methods))
    smry[["methods"]] <- smry[["methods"]][methods]
  
  class(smry) <- c("boolsum", class(smry))
  smry
}


##' Print summary of boolean model as described in \code{boolsum} object.
##'
##' This function prints a summary of a boolean model as described in a
##' \code{boolsum}, itself the result of calling \code{summary} on an object of
##' class \code{boolean}.
##'
##' @method print boolsum
##' @param x object of class \code{boolsum}.
##' @param ... Additional parameters to be passed.
##' @return Prints summary information about the boolean model fit.
##' @author Jason W. Morgan (\email{morgan.746@@osu.edu})
##' @export print.boolsum
print.boolsum <- function(x, ...) {
  cat("\n", paste(rep("=", 74), collapse=""), "\n",
      "Boolean model estimate\n\n", sep="")

  cat(paste(format(c("Observations", "Parameters")),
            paste(format(c(x$est.smry[1], x$est.smry[2])),
                  "\n")), sep="")
  
  for (i in names(x$methods)) {
    X <- x$methods[[i]]
    cat("\n", paste(rep("-", 30), collapse=""), "\n",
        "Estimation method: ", i, "\n",
        paste(rep("-", 30), collapse=""), "\n", sep="")

    print(format(round(X$coef, 3), width=8, scientific=FALSE), quote=FALSE,
          right=TRUE)
    cat("\n")

    if (i == "mcmc")
      print(format(X$method.smry), quote=FALSE)
    else
      print(format(X$method.smry, nsmall=3), quote=FALSE)
  }
  
  cat(paste(rep("=", 74), collapse=""), "\n\n")
  misc_msg(1)  
}


##' Default print method for boolprof objects.
##'
##' This is the default print object for \code{boolprof} objects.
##' 
##' @method print boolprof
##' @param x object of the \code{boolprof-class}. 
##' @param ... Additional arguments that are passed to the \code{lattice} plot
##' object.
##' @return Plots the estimated log-likelihood profiles.
##' @author Jason W. Morgan (\email{morgan.746@@osu.edu})
##' @export print.boolprof
print.boolprof <- function(x, ...) {
  plot(x$default.plot, ...)
}


##' Default print method for boolprob objects.
##'
##' This is the default print object for \code{boolprob} objects.
##' 
##' @method print boolprob
##' @param x object of the \code{boolprob-class}. 
##' @param ... Additional arguments that are passed to the \code{lattice} plot
##' object.
##' @return Plots the estimated predicted probabilities.  
##' @author Jason W. Morgan (\email{morgan.746@@osu.edu})
##' @export print.boolprob
print.boolprob <- function(x, ...) {
  plot(x$default.plot, ...)
}


##' Profile Likelihood Plot
##'
##' Plot the profile likelihoods for selected covariates. Adjust default plot
##' produced by \code{\link{boolprof}} to fit the user's preferences.
##'
##' @method plot boolprof
##' @param x boolprof object. 
##' @param y character or numeric vector specifying the covariates for which
##' the profiled likelihood should be plotted. 
##' @param ... additional arguments to pass to  \code{\link{update.trellis}}.
##' See \code{\link{xyplot}} for details.
##' @return Plots the profile likelihood for the selected covariates.
##' @author Jason W. Morgan <morgan.746@@osu.edu>
##' @export plot.boolprof
plot.boolprof <- function(x, y=NULL, ...) {
  if (is.null(y)) {
    print(x, ...)
  }
  else {
    if (is.character(y))
      y <- match(y, x$coef.labels)
    
    update(x$default.plot[y], ...)
  }
}


##' Profile Predicted Probability Plot
##'
##' Plot predicted probabilities for selected covariates. Adjust default plot
##' produced by \code{\link{boolprob}} to fit the user's preferences.
##'
##' @method plot boolprob
##' @param x boolprob object. 
##' @param y character or numeric vector specifying the covariates for which
##' the predicted probabilities should be plotted. 
##' @param ... additional arguments to pass to  \code{\link{update.trellis}}.
##' See \code{\link{xyplot}} for details.
##' @return Plots the predicted probabilities for the selected covariates.
##' @author Jason W. Morgan (\email{morgan.746@@osu.edu})
##' @export plot.boolprob
plot.boolprob <- function(x, y=NULL, ...) {
  if (is.null(y)) {
    print(x, ...)
  }
  else {
    if (is.character(y))
      y <- match(y, x$coef.labels)
    
    update(x$default.plot[y], ...)
  }
}


##' Default print for boolboot objects.
##'
##' Default print for boolboot objects.
##' @method print boolboot
##' @title Print Bootstrap Results for Boolean Object 
##' @param x boolboot object.
##' @param ... Additional parameters to pass on.
##' @return Print results in abbreviated summary.
##' @author Jason W. Morgan (\email{morgan.746@@osu.edu})
##' @export print.boolboot
print.boolboot <- function(x, ...) {
  beta <- x$model.boot$coef
  n <- x$model.boot$samples
  e <- n - nrow(beta)
  
  ## Create summary
  b <- apply(beta, 2, function(x) {
    c(mean(x, na.rm=TRUE), median(x, na.rm=TRUE), sd(x, na.rm=TRUE),
      quantile(x, probs=c(0.025, 0.975), na.rm=TRUE))
  })

  est <- round(t(b), 3)
  colnames(est) <- c("Mean", "Median", "SD", "2.5%", "97.5%")
  rownames(est) <- x$coef.labels

  ## Print results
  cat("\n")
  cat(paste(rep("=", 74), collapse=""), "\n")
  cat("Bootstrap results from boolean model estimate\n\n")
  print(est)
  cat("\nNumber of bootstrap samples:", n, "\n")
  ##cat("Number of failed convergences:", e, "\n")
  if (n < 1000) misc_msg(61)
  ## if (e / n >= 0.1) misc_msg(62)
  ## cat(paste(rep("=", 74), collapse=""), "\n\n")
}

##' Summary function for boolboot objects. 
##'
##' Summary function for boolboot objects.
##' @method summary boolboot
##' @title Print Bootstrap Results for Boolean Object 
##' @param object object of class boolboot.
##' @param ... additional parameters to pass on.
##' @return Prints a summary.
##' @author Jason W. Morgan (\email{morgan.746@@osu.edu})
##' @export summary.boolboot
summary.boolboot <- function(object, ...) {
  print.boolboot(object, ...)
}
