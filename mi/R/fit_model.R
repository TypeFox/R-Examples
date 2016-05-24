# Part of the mi package for multiple imputation of missing data
# Copyright (C) 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015 Trustees of Columbia University
# 
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.


## these fit a regression and return the model

# note, helper functions are good because they are checked more rigorously by R CMD check

setMethod("fit_model", signature(y = "missing_variable", data = "missing_data.frame"), def =
  function(y, data, s, warn, ...) {
    stop("This method should not have been called. You need to define the relevant fit_model() S4 method")
  })

setMethod("fit_model", signature(y = "irrelevant", data = "missing_data.frame"), def = 
  function(y, data, ...) {
    stop("'fit_model' should not have been called on an 'irrelevant' variable")
  })

setMethod("fit_model", signature(y = "binary", data = "missing_data.frame"), def =
  function(y, data, s, warn, X = NULL, ...) {
    if(is.null(X)) {
      to_drop <- data@index[[y@variable_name]]
      if(length(to_drop)) X <- data@X[,-to_drop]
      else X <- data@X[,]
      if(is(data, "experiment_missing_data.frame")) {
        treatment <- names(which(data@concept == "treatment"))
        if(data@concept[y@variable_name] == "outcome") {
          X <- cbind(X, interaction = X * data@variables[[treatment]]@data)
        }
      }
    }
    
    if(s > 1)       start <- y@parameters[s-1,]
    else if(s < -1) start <- y@parameters[1,]
    else start <- NULL
    start <- NULL
    weights <- if(length(data@weights) == 1) data@weights[[1]] else data@weights[[y@variable_name]]
    CONTROL <- list(epsilon = max(1e-8, exp(-abs(s))), maxit = 25, trace = FALSE)
    priors <- data@priors[[y@variable_name]]
    out <- bayesglm.fit(X, y@data - 1L, weights = weights,
                                   prior.mean = priors[[1]], prior.scale = priors[[2]], prior.df = priors[[3]],
                                   prior.mean.for.intercept = priors[[4]], prior.scale.for.intercept = priors[[5]],
                                   prior.df.for.intercept = priors[[6]],
                                   start = start, family = y@family, Warning = FALSE, control = CONTROL)
    if(warn && !out$converged) {
      warning(paste("bayesglm() did not converge for variable", y@variable_name, "on iteration", abs(s)))
    }
    if(any(abs(coef(out)) > 100)) {
      warning(paste(y@variable_name, ": separation on iteration", abs(s)))
    }
    out$x <- X
    class(out) <- c("bayesglm", "glm", "lm")
    return(out)
  })

.fit_MNL <-
  function(y, X, weights) {
    model<-nnet::multinom(y@data ~ X -1, weights = weights, Hess = y@imputation_method == "ppd", 
                           model = TRUE, trace = FALSE, MaxNWts = 10000)
    return(model)
  }

.fit_RNL <-
  function(y, X, weights, CONTROL) {
    if (y@use_NA==TRUE) values <- c(-.Machine$integer.max, 1:length(y@levels))  
    else values <- 1:length(y@levels)
    out <- sapply(values, simplify = FALSE, FUN = function(l) {
      model <- bayesglm.fit(X, y@data == l, weights = weights, 
                            family = y@family, control=CONTROL)
      model$x <- X # bayesglm.fit() by default does not retain the model matrix it uses 
      class(model) <- c("bayesglm", "glm", "lm")  
      return(model)
    })
    class(out) <- "RNL"
    return(out)
  }

setMethod("fit_model", signature(y = "unordered-categorical", data = "missing_data.frame"), def =
  function(y, data, warn, s, ...) {
    to_drop <- data@index[[y@variable_name]]
    if (y@use_NA) {
      y@data[y@which_drawn] <- -.Machine$integer.max # make NAs the smallest possible signed integer 
    }
	if(length(to_drop)) X <- data@X[,-to_drop]
    else X <- data@X[,]
    if(is(data, "experiment_missing_data.frame")) {
      treatment <- names(which(data@concept == "treatment"))
      if(data@concept[y@variable_name] == "outcome") {
        X <- cbind(X, interaction = X * data@variables[[treatment]]@data)
      }
    }
    weights <- if(length(data@weights) == 1) data@weights[[1]] else data@weights[[y@variable_name]]
    if(y@estimator == "MNL") {
      out <- .fit_MNL(y, X, weights)
      data@X
    }
  	else if(y@estimator == "RNL"){
      CONTROL <- list(epsilon = max(1e-8, exp(-abs(s))), maxit = 25, trace = FALSE)
      out <- .fit_RNL(y, X, weights, CONTROL)
      data@X
    }
    else stop("estimator not recognized")
    return(out)
  })

.clogit <- # similar to the survival::clogit function
  function(formula, data, n, method, weights, subset, 
           x = TRUE, na.action = "na.exclude") {
    coxcall <- match.call()
    coxcall[[1]] <- as.name("coxph")
    newformula <- formula
    newformula[[2]] <- substitute(survival::Surv(rep(1, nn), case), 
                                  list(case = formula[[2]], nn = n))
    environment(newformula) <- environment(formula)
    coxcall$formula <- newformula
    coxcall$n <- NULL
    coxcall <- eval(coxcall, sys.frame(sys.parent()))
    coxcall$userCall <- sys.call()
    class(coxcall) <- c("clogit", "coxph")
    coxcall
  }

setMethod("fit_model", signature(y = "grouped-binary", data = "missing_data.frame"), def =
  function(y, data, s, warn) {
    # see http://www.stata.com/support/faqs/stat/clogitcl.html for a good explanation of this model
    to_drop <- data@index[[y@variable_name]]
    X <- data@X[,-to_drop]
    weights <- if(length(data@weights) == 1) data@weights[[1]] else data@weights[[y@variable_name]]
    groups <- sapply(y@strata, FUN = function(x) complete(data@variables[[x]], m = 0L), simplify = FALSE)  
    out <- .clogit(y@data ~ X + strata(groups), method = "breslow", weights = weights, n = nrow(X)) 
    out$x <- X
    return(out)
  })

setMethod("fit_model", signature(y = "ordered-categorical", data = "missing_data.frame"), def =
  function(y, data, s, warn, X = NULL, ...) {
    if(is.null(X)) {
      to_drop <- data@index[[y@variable_name]]
      if(length(to_drop)) X <- data@X[,-to_drop]
      else X <- data@X[,]
      if(is(data, "experiment_missing_data.frame")) {
        treatment <- names(which(data@concept == "treatment"))
        if(data@concept[y@variable_name] == "outcome") {
          X <- cbind(X, interaction = X * data@variables[[treatment]]@data)
        }
      }
      X <- X[,-1]
    }
    
    method <- if(y@family$link == "logit") "logistic" else y@family$link
    start <- NULL
    start <- c(rep(0, ncol(X)), qlogis(cumsum(table(y@data)) / nrow(X)))
    start <- start[-length(start)]
    weights <- if(length(data@weights) == 1) data@weights[[1]] else data@weights[[y@variable_name]]
    CONTROL <- list(reltol = max(1e-8, exp(-abs(s))))
    priors <- data@priors[[y@variable_name]]
    out <- bayespolr(as.ordered(y@data) ~ X, weights = weights, method = method,
                                prior.mean = priors[[1]], prior.scale = priors[[2]], prior.df = priors[[3]],
                                prior.counts.for.bins = priors[[4]], control = list(reltol = max(1e-8, exp(-abs(s)))), ...)
    
    if(warn && out$convergence != 0) {
      warning(paste("bayespolr() did not converge for variable", y@variable_name, "on iteration", abs(s)))
    }
    out$x <- X
    return(out)
  })

setMethod("fit_model", signature(y = "interval", data = "missing_data.frame"), def =
  function(y, data, s, warn, ...) {
    stop("FIXME: write this method")
  })

## helper function
.fit_continuous <-
  function(y, data, s, warn, X, subset = 1:nrow(X)) {
    weights <- if(length(data@weights) == 1) data@weights[[1]] else data@weights[[y@variable_name]]
    if(!is.null(weights)) weights <- weights[subset]
    
    if(s > 1)       start <- y@parameters[s-1,]
    else if(s < -1) start <- y@parameters[1,]
    else start <- NULL
    start <- NULL
    mark <- c(TRUE, apply(X[subset,-1, drop = FALSE], 2, FUN = function(x) length(unique(x)) > 1))
    if(!all(mark)) {
      if(abs(s) == 1) {
        stop(paste(y@variable_name, ": imputed values on iteration 0 randomly inadmissible; try mi() again with different seed"))
      }
      X <- X[,mark]
      if(!is.null(start)) start <- start[mark]
    }
    CONTROL <- list(epsilon = max(1e-8, exp(-abs(s))), maxit = 25, trace = FALSE)
    priors <- data@priors[[y@variable_name]]
    out <- bayesglm.fit(X[subset,], y@data[subset], weights = weights, start = start, family = y@family, 
                                   prior.mean = priors[[1]], prior.scale = priors[[2]], prior.df = priors[[3]],
                                   prior.mean.for.intercept = priors[[4]], prior.scale.for.intercept = priors[[5]],
                                   prior.df.for.intercept = priors[[6]], Warning = FALSE, control = CONTROL)
    
    if(warn && !out$converged) {
      warning(paste("bayesglm() did not converge for variable", y@variable_name, "on iteration", abs(s)))
    }
    out$x <- X
    class(out) <- c("bayesglm", "glm", "lm")
    return(out)
  }

setMethod("fit_model", signature(y = "continuous", data = "missing_data.frame"), def =
  function(y, data, s, warn, ...) {
    to_drop <- data@index[[y@variable_name]]
    if(length(to_drop)) X <- data@X[,-to_drop]
    else X <- data@X[,]
    if(is(data, "experiment_missing_data.frame")) {
      treatment <- names(which(data@concept == "treatment"))
      if(data@concept[y@variable_name] == "outcome") {
        X <- cbind(X, interaction = X * data@variables[[treatment]]@data)
      }
    }
    return(.fit_continuous(y, data, s, warn, X))
  })

# setMethod("fit_model", signature(y = "truncated-continuous", data = "missing_data.frame"), def =
# function(y, data, s, warn, ...) {
#   stop("FIXME: write this method using library(survival)")
# })
# 
# setMethod("fit_model", signature(y = "censored-continuous", data = "missing_data.frame"), def =
# function(y, data, s, warn, ...) {
#   stop("FIXME: mi does not do censored-continuous variables yet")
#   to_drop <- data@index[[y@variable_name]]
#   X <- cbind(y@raw_data, data@X[,-to_drop])
#   if(is(data, "experiment_missing_data.frame")) {
#     treatment <- names(which(data@concept == "treatment"))
#     if(data@concept[y@variable_name] == "outcome") {
#       X <- cbind(X, interaction = X * data@variables[[treatment]]@data)
#     }
#   }  
# })

setMethod("fit_model", signature(y = "semi-continuous", data = "missing_data.frame"), def =
  function(y, data, s, warn, ...) {
    stop("the semi-continuous fit_model() method should not have been called")
  })

setMethod("fit_model", signature(y = "nonnegative-continuous", data = "missing_data.frame"), def =
  function(y, data, s, warn, ...) {  
    to_drop <- data@index[[y@variable_name]]
    if(length(to_drop)) X <- data@X[,-to_drop]
    else X <- data@X[,]
    model <- fit_model(y@indicator, data, s, warn, X)
    if(abs(s) > 1) subset <- complete(y@indicator, m = 0L, to_factor = TRUE) == 0
    else subset <- 1:nrow(X)
    return(.fit_continuous(y = y, data = data, s = s, warn = warn, X = X, subset = subset))
  })

setMethod("fit_model", signature(y = "SC_proportion", data = "missing_data.frame"), def =
  function(y, data, s, warn, ...) {
    to_drop <- data@index[[y@variable_name]]
    if(length(to_drop)) X <- data@X[,-to_drop]
    else X <- data@X[,]
    model <- fit_model(y@indicator, data, s, warn, X)
    if(abs(s) > 1) subset <- complete(y@indicator, m = 0L, to_factor = TRUE) == 0
    else subset <- 1:nrow(X)
    return(.fit_proportion(y = y, data = data, s = s, warn = warn, X = X, subset = subset))
  })

## helper function
.fit_proportion <-
  function(y, data, s, warn, X, subset = 1:nrow(X)) {
    weights <- if(length(data@weights) == 1) data@weights[[1]] else data@weights[[y@variable_name]]
    if(!is.null(weights)) weights <- weights[subset]
    
    if(s > 1)       start <- y@parameters[s-1,]
    else if(s < -1) start <- y@parameters[1,]
    else start <- NULL
    start <- NULL
    mark <- c(TRUE, apply(X[subset,-1, drop = FALSE], 2, FUN = function(x) length(unique(x)) > 1))
    if(!all(mark)) {
      if(abs(s) == 1) {
        stop(paste(y@variable_name, ": imputed values on iteration 0 randomly inadmissible; try mi() again with a different seed"))
      }
      X <- X[,mark]
      if(!is.null(start)) start <- start[c(mark, TRUE)]
    }
    out <- betareg::betareg.fit(X[subset,], y@data[subset], weights = if(!is.null(weights)) weights[subset],
                                      link = y@family$link, link.phi = y@link.phi, 
                                      control = betareg::betareg.control(reltol = 1e-8, start = start, fsmaxit = 0))
    if(warn && !out$converged) {
      warning(paste("betareg() did not converge for variable", y@variable_name, "on iteration", abs(s)))
    }
    out$x <- X
    class(out) <- c("betareg")
    return(out)
  }

setMethod("fit_model", signature(y = "proportion", data = "missing_data.frame"), def =
  function(y, data, s, warn, ...) {
    to_drop <- data@index[[y@variable_name]]
    if(length(to_drop)) X <- data@X[,-to_drop]
    else X <- data@X[,]
    if(y@family$family == "gaussian") out <- .fit_continuous(y, data, s, warn, X)
    else out <- .fit_proportion(y, data, s, warn, X)
    return(out)
  })

setMethod("fit_model", signature(y = "count", data = "missing_data.frame"), def =
  function(y, data, s, warn, ...) {
    to_drop <- data@index[[y@variable_name]]
    if(length(to_drop)) X <- data@X[,-to_drop]
    else X <- data@X[,]
    return(.fit_continuous(y, data, s, warn, X)) # even though counts are not continuous
  })

## experiments
setMethod("fit_model", signature(y = "missing_variable", data = "experiment_missing_data.frame"), def =
  function(y, data, ...) {
    stop("you need to write a specific fit_model() method for the", class(y), "class")
  })

setMethod("fit_model", signature(y = "continuous", data = "experiment_missing_data.frame"), def =
  function(y, data, s, warn, ...) {
    to_drop <- data@index[[y@variable_name]]
    ## For each case, make an X matrix based on the giant matrix in data@X
    if(data@case == "outcomes") { # missingness on outcome(s) only
      if(length(to_drop)) X <- data@X[,-to_drop]
      else X <- data@X[,]
      treatment_name <- names(data@concept[data@concept == "treatment"])
      X <- cbind(X, interaction = X[,!(colnames(X) %in% c("(Intercept)", treatment_name))] * X[,treatment_name])
    }
    else if(data@case == "covariates") { # missingness on covariate(s) only
      to_drop <- c(to_drop, which(data@concept == "treatment"))
      if(length(to_drop)) X <- data@X[,-to_drop]
      else X <- data@X[,]
    }
    else { # missing on both outcome(s) and covariate(s)
      if(data@concept[y@variable_name] == "covariate") {
        to_drop <- c(to_drop, which(data@concept == "treatment"))
      }
      if(length(to_drop)) X <- data@X[,-to_drop]
      else X <- data@X[,]
    }
    return(mi::.fit_continuous(y, data, s, warn, X))
  })

## here y indicates which variable to model
setMethod("fit_model", signature(y = "character", data = "mi"), def = 
  function(y, data, m = length(data@data), ...) {
    s <- sum(data@total_iters) + 1
    if(length(m) == 1) {
      models <- vector("list", m)
      for(i in 1:m) {
        model <- data@data[[i]]@variables[[y]]@model
        if(is.null(model)) {
          model <- fit_model(y = data@data[[i]]@variables[[y]], data = data@data[[i]], s = s, warn = TRUE, ...)
          if(!isS4(model)) model$x <- model$X <- model$y <- NULL
        }
        models[[i]] <- model
      }
    }
    else {
      models <- vector("list", length(m))
      models <- for(i in 1:length(m)) {
        if(is.null(data@data[[i]]@variables[[y]]@model)) {
          models[[m[i]]] <- fit_model(y = data@data[[m[i]]]@variables[[y]], data = data@data[[m[i]]], s = s, warn = TRUE, ...)
        }
        else models[[i]] <- data@data[[i]]@variables[[y]]@model
      }
    }
    return(models)
  })

## fit all variables with missingness
setMethod("fit_model", signature(y = "missing", data = "mi"), def = 
  function(data, m = length(data@data)) {
    varnames <- names(data@data[[1]]@variables)
    exclude <- data@data[[1]]@no_missing | 
               sapply(data@data[[1]], FUN = function(y) is(y, "irrelevant"))
    models <- sapply(varnames, simplify = FALSE, FUN = function(v) {
                if(v %in% exclude) paste(v, "not modeled") ## maybe just skip these? 
                else fit_model(y = v, data = data, m = m)
                })
    return(models)
  })

## fit all elements of a mdf_list
setMethod("fit_model", signature(y = "missing", data = "mdf_list"), def =
  function(data, s = -1, verbose = FALSE, warn = FALSE, ...) {
    out <- lapply(data, fit_model, s = s, verbose = verbose, warn = warn, ...)
    class(out) <- "mdf_list"
    return(out)
  })

.fit_model_y <-
  function(y, data, s, verbose, warn, ...) {
    if(s != 0 && y@imputation_method != "mcar") {
      if(is(y, "semi-continuous")) {
        to_drop <- data@index[[y@variable_name]]
        if(length(to_drop)) X <- data@X[,-to_drop]
        else X <- data@X[,]
        model <- fit_model(y = y@indicator, data = data, s = s, warn = warn, X = X)
        indicator <-    mi(y = y@indicator, model = model, s = ifelse(s < 0, 1L, s))
        if(s > 1) indicator@parameters[s,] <- get_parameters(model)
        else if(abs(s) == 1) {
          parameters <- get_parameters(model)
          rows <- if(s == 1) nrow(indicator@parameters) else 1
          if(ncol(indicator@parameters) == 0) {
            temp <- matrix(NA_real_, nrow = rows, ncol = length(parameters))
          }
          temp[1,] <- parameters
          indicator@parameters <- temp
        }
        else indicator@parameters[1,] <- get_parameters(model)
        y@indicator <- indicator
      }
      model <- fit_model(y = y, data = data, s = s, warn = warn)
      y <- mi(y = y, model = model, s = ifelse(s < 0, 1L, s))
    }
    else y <- mi(y = y)
    if(y@imputation_method == "mcar") {
      # do nothing
    }
    else if(s > 1) {
      parameters <- get_parameters(model)
      if(length(parameters) != ncol(y@parameters)) parameters <- y@parameters[s-1,] # scary
      y@parameters[s,] <- parameters
    }
    else if(abs(s) == 1) {
      parameters <- get_parameters(model)
      rows <- if(s == 1) nrow(y@parameters) else 1
      if(ncol(y@parameters) == 0) {
        temp <- matrix(NA_real_, nrow = rows, ncol = length(parameters))
      }
      temp[1,] <- parameters
      y@parameters <- temp
    }
    else if(s != 0) {
      parameters <- get_parameters(model)
      if(length(parameters) == ncol(y@parameters)) y@parameters[s,] <- parameters
    }
    
    return(y)
  }

.update_X <-
  function(y, data) {
    which_drawn <- y@which_drawn
    varname <- y@variable_name
    if(is(y, "categorical")) {
      dummies <- .cat2dummies(y)[which_drawn,,drop = FALSE]
      data@X[ which_drawn, data@index[[varname]][1:NCOL(dummies)]] <- dummies
    }
    else if(is(y, "semi-continuous")) {
      mark <- data@index[[varname]]
      data@X[ which_drawn, mark[1]  ] <- y@data[which_drawn]
      dummies <- .cat2dummies(y@indicator)
      data@X[ which_drawn, mark[1 + 1:NCOL(dummies)] ] <- dummies[which_drawn,,drop = FALSE]
    }
    else if(is(y, "censored_continuous")) {
      temp <- y@data[which_drawn]
      if(y@n_lower) temp <- cbind(temp, lower = y@lower_indicator@data[y@which_drawn])
      if(y@n_upper) temp <- cbind(temp, upper = y@upper_indicator@data[y@which_drawn])
      data@X[ which_drawn, data@index[[varname]][1:NCOL(temp)]] <- temp
      data@X[ y@which_censored, data@index[[varname]][1] ] <- y@data[y@which_censored]
    }
    else if(is(y, "truncated_continuous")) {
      temp <- y@data[which_drawn]
      if(y@n_lower) temp <- cbind(temp, lower = y@lower_indicator@data[y@which_drawn])
      if(y@n_upper) temp <- cbind(temp, upper = y@upper_indicator@data[y@which_drawn])
      data@X[ which_drawn, data@index[[varname]][1:NCOL(temp)]] <- temp
      data@X[ y@which_truncated, data@index[[varname]][1] ] <- y@data[y@which_truncated]
    }
    else data@X[ which_drawn, data@index[[varname]][1] ] <- y@data[which_drawn]
    return(data)    
  }

.fit_model_mdf <-
  function(data, s, verbose, warn, ...) {
    if(verbose) {
      txt <- paste("Iteration:", abs(s))
      if(isatty(stdout()) && !(any(search() == "package:gWidgets"))) cat("\n", txt)
      else cat("<br>", txt, file = file.path(data@workpath, "mi.html"), append = TRUE)
    }
    on.exit(print("the problematic variable is"))
    on.exit(show(y), add = TRUE)
    for(jj in sample(1:ncol(data), ncol(data), replace = FALSE)) {
      y <- data@variables[[jj]]
      if(y@all_obs) next
      if(is(y, "irrelevant")) next
      y <- .fit_model_y(y, data, s, verbose, warn, ...)
      data <- .update_X(y, data)
      data@variables[[jj]] <- y
      if(verbose) {
        txt <- "."
        if(isatty(stdout()) && !(any(search() == "package:gWidgets"))) cat(txt)
        else cat(txt, file = file.path(data@workpath, "mi.html"), append = TRUE)
      }
    }
    if(verbose) {
      if(isatty(stdout()) && !(any(search() == "package:gWidgets"))) cat(" ")
      else cat(" <br/>", file = file.path(data@workpath, "mi.html"), append = TRUE)
    }
    if(.MI_DEBUG) sapply(data@variables, validObject, complete = TRUE)
    on.exit()
    return(data)
  }

## unlike the above methods, these return a (modified) missing_data.frame
setMethod("fit_model", signature(y = "missing", data = "missing_data.frame"), def =
  function(data, s = -1, verbose = FALSE, warn = FALSE, ...) {
    return(.fit_model_mdf(data = data, s = s, verbose = verbose, warn = warn, ...))
  })

setMethod("fit_model", signature(y = "missing", data = "allcategorical_missing_data.frame"), def =
  function(data, s = -1, verbose = FALSE, warn = FALSE, ...) {
    if(verbose) {
      txt <- paste("Iteration:", abs(s))
      if(isatty(stdout()) && !(any(search() == "package:gWidgets"))) cat("\n", txt)
      else cat("<br>", txt, file = file.path(data@workpath, "mi.html"), append = TRUE)
    }
    
    Hstar <- data@Hstar
    
    if(abs(s) == 0) { # starting iteration
      V_h <- c(runif(Hstar - 1), 1)
      c_prod <- cumprod(1 - V_h)
      data@parameters$pi <- V_h * c(1, c_prod[-Hstar])
      data@variables <- lapply(data@variables, FUN = function(y) {
        if(is(y, "irrelevant")) return(y)
        if(y@all_obs) return(y)
        return(mi(y)) # bootstrapping
      })
      data@X <- do.call(cbind, args = lapply(data@variables, FUN = function(y) {
        if(is(y, "irrelevant")) return(NULL)
        else return(y@data)
      }))
      phi <- lapply(1:Hstar, FUN = function(h) {
        lapply(data@variables, FUN = function(y) {
        if(is(y, "irrelevant")) return(NULL)
        return(c(tabulate(y@data, nbins = length(y@levels)) / y@n_total))
      })})
      data@parameters$phi <- phi
      data@parameters$alpha <- 1
      
      cols <- Hstar
      rows <- nrow(data@latents@imputations)
      if(ncol(data@latents@parameters) == 0) {
        temp <- matrix(NA_real_, nrow = rows, ncol = cols)
      }
      data@latents@parameters <- temp
      
      return(data)
    }
      
    # S1: Update latent class membership
    pi  <- data@parameters$pi
    phi <- data@parameters$phi
    probs <- sapply(1:Hstar, FUN = function(h) {
           phi_h <- phi[[h]]
           numerators <- rep(1, nrow(data))
           for(j in 1:ncol(data)) {
             y <- data@variables[[j]]
             if(is(y, "irrelevant")) next
             phi_hj <- phi_h[[y@variable_name]]
             numerators <- numerators * 
                           data@X[,y@variable_name] * phi_hj[data@X[,y@variable_name]]
           }
           numerators <- numerators * pi[h]
           return(numerators)
    })
    z <- apply(probs, 1, FUN = function(prob) {
                 which(rmultinom(1,1,prob) == 1) # rmultinom normalizes internally
              })
    data@latents@data[] <- z
    data@latents@imputations[s,] <- z
    data@latents@parameters[s,] <- pi
    
    # S2: Update V_h
    n_h <- c(tabulate(z, nbins = Hstar))
    V_h <- sapply( 1:(Hstar - 1), FUN = function(h) {
                   a <- 1 + n_h[h]
                   b <- data@parameters$alpha + sum(n_h[-c(1:h)])
                   if(b == 0) return(1)
                   rbeta(1, a, b) 
                  })
    V_h <- c(V_h, 1)
    c_prod <- cumprod(1 - V_h)
    data@parameters$pi <- V_h * c(1, c_prod[-Hstar])
    
    # S3: Update choice probabilities
    phi <- lapply(1:Hstar, FUN = function(h) lapply(data@variables, FUN = function(y) {
      if(is(y, "irrelevant")) return(NULL)
      mark <- z == h
      tab <- tabulate(y@data[mark], nbins = length(y@levels))
      return(.rdirichlet(1, data@priors$a[y@variable_name] + c(tab)))
    }))
    data@parameters$phi <- phi
    
    # S4: Update alpha
    alpha <- rgamma(1, data@priors$a_alpha + Hstar - 1, 
                    data@priors$b_alpha - log(pi[Hstar]))
    data@parameters$alpha <- alpha
    
    # S5: Impute
    data@variables <- lapply(data@variables, FUN = function(y) {
      if(is(y, "irrelevant")) return(y)
      if(y@all_obs) return(y)
      if(verbose) {
        txt <- "."
        if(isatty(stdout()) && !(any(search() == "package:gWidgets"))) cat(txt)
        else cat(txt, file = file.path(data@workpath, "mi.html"), append = TRUE)
      }      
      classes <- z[y@which_drawn]
      uc <- unique(classes)
      Pr <- t(sapply(uc, FUN = function(c) phi[[c]][[y@variable_name]]))
      rownames(Pr) <- uc
      y <- mi(y, Pr[as.character(classes),,drop=FALSE])
    })
    data@X <- do.call(cbind, args = lapply(data@variables, FUN = function(y) {
      if(is(y, "irrelevant")) return(NULL)
      else return(y@data)
    }))

    return(data)
})

.fit_model_Sophie <-
  function(y, data, s = -1, verbose = FALSE, warn = FALSE, ...) {
    classes <- data@latents@data
    uc <- unique(classes)
    Pr <- t(sapply(uc, FUN = function(c) data@parameters$phi[[c]][[y@variable_name]]))
    rownames(Pr) <- uc
#     Pr <- Pr[as.character(classes),,drop=FALSE]
    Pr <- Pr / rowSums(Pr)
    return(list(fitted = Pr))    
  }

setMethod("fit_model", signature(y = "unordered-categorical", 
                                 data = "allcategorical_missing_data.frame"), def =
  function(y, data, s = -1, verbose = FALSE, warn = FALSE, ...) {
    return(.fit_model_Sophie(y, data, s, verbose, warn, ...))
  })

setMethod("fit_model", signature(y = "ordered-categorical", 
                                 data = "allcategorical_missing_data.frame"), def =
  function(y, data, s = -1, verbose = FALSE, warn = FALSE, ...) {
    return(.fit_model_Sophie(y, data, s, verbose, warn, ...))
})

setMethod("fit_model", signature(y = "binary", 
                                 data = "allcategorical_missing_data.frame"), def =
  function(y, data, s = -1, verbose = FALSE, warn = FALSE, ...) {
    return(.fit_model_Sophie(y, data, s, verbose, warn, ...))
})

setMethod("fit_model", signature(y = "missing_data.frame", data = "missing_data.frame"), def =
  function(y, data, s = -1, verbose = FALSE, warn = FALSE, ...) {
    
    if(verbose) {
      txt <- paste("Iteration:", abs(s))
      if(isatty(stdout()) && !(any(search() == "package:gWidgets"))) cat("\n", txt)
      else cat("<br>", txt, file = file.path(data@workpath, "mi.html"), append = TRUE)
    }
    for(jj in sample(1:ncol(y), ncol(y), replace = FALSE)) {
      z <- y@variables[[jj]]
      if(z@all_obs) next
      if(is(z, "irrelevant")) next
      y@variables[[jj]] <- .fit_model_y(z, data, s, verbose, warn, ...)
      if(verbose) {
        txt <- "."
        if(isatty(stdout()) && !(any(search() == "package:gWidgets"))) cat(txt)
        else cat(txt, file = file.path(data@workpath, "mi.html"), append = TRUE)
      }
    }
    if(verbose) {
      if(isatty(stdout()) && !(any(search() == "package:gWidgets"))) cat(" ")
      else cat(" <br/>", file = file.path(data@workpath, "mi.html"), append = TRUE)
    }
    if(.MI_DEBUG) sapply(data@variables, validObject, complete = TRUE)
    return(y)
    
  })

setMethod("fit_model", signature(y = "missing", data = "multilevel_missing_data.frame"), def =
  function(data, s = -1, verbose = FALSE, warn = FALSE, ...) {
    data@mdf_list <- fit_model(data = data@mdf_list, s = s, verbose = verbose, warn = warn, ...)
    if(s == 0) return(data)
    
    ## FIXME: Implement 3+ levels recursively
    
    # update group means
    means <- sapply(data@mdf_list, FUN = function(x) colMeans(x@X[,-1]))
    if(is.list(means)) {
      
    }
    else means <- t(means)
    mark <- 0L ## FIXME
    data@X[,mark] <- means
    
    # impute the group level variables if necessary
    data <- .fit_model_mdf(data = data, s = s, verbose = verbose, warn = warn, ...)
    
    # model the individual level estimates
    for(i in seq_along(ncol(data))) {
      if(is(data@variables[[i]], "irrelevant")) next
      #     if(data@no_missing[i]) next
      mark <- if(s < 0) 1 else s
      fish <- sapply(data@mdf_list, FUN = function(d) d@variables[[i]]@parameters[mark,])
      if(is.list(fish)) {
        ## FIXME: may be a list
      }
      else fish <- t(fish)
      
      for(j in seq_along(ncol(fish))) {
        model <- bayesglm.fit(data@X, y = fish[,j]) # group-level regression
        class(model) <- c("bayesglm", "glm", "lm")
        params <- arm::sim(model, 1)
        beta <- params@coef
        sigma <- params@sigma
        yhats <- rnorm(nrow(X), data@X %*% beta, sd = sigma)
        # change the priors for each element of the mdf_list accordingly
        for(k in seq_along(data@mdf_list)) {
          data@mdf_list[[k]]$mean[colnames(data)[j]] <- yhats[k]
          data@mdf_list[[k]]$sd[colnames(data)[j]] <- sigma
        }
      }
    }
    return(data)
  })
