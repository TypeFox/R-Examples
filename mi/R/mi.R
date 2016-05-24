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

.prune_missing_variable <-
  function(y, s) {
    if(!is(y, "missing_variable")) stop("'y' must inherit from the 'missing_variable' class")
    if(!y@all_obs) {
      y@parameters <- y@parameters[1:s,,drop = FALSE]
      y@imputations <- y@imputations[1:s,,drop = FALSE]
    }
    return(y)
  }

.MPinverse <- function(eta, tol = sqrt(.Machine$double.eps)) {
  cov_eta <- cov(eta)
  ev <- eigen(cov_eta, TRUE)
  ev$values <- ifelse(ev$values > tol, 1/ev$values, 0)
  Sigma_inv <- crossprod(sqrt(ev$values)*(t(ev$vectors)))
  return(Sigma_inv)
}

.mi <- function(i, y, verbose, s_start, s_end, ProcStart, max.minutes, parallel, save_models) {
  mdf <- y
  if(verbose) message("Chain ", i, "\n")
  for(s in s_start:s_end) {
    if(verbose) message("Chain ", i, " Iteration ", s, "\n")
    mdf <- fit_model(data = mdf, s = s, verbose = FALSE, warn = s == s_end)
    if(s > 0) {
      pars <- unlist(sapply(mdf@variables, FUN = function(y) {
        if(is(y, "irrelevant")) return(NA_real_)
        else if(y@all_obs) return(NA_real_)
        else return(y@parameters[s,,drop=TRUE])
      }))
      pars <- t(pars[!is.na(pars)])
      fp <- file.path(mdf@workpath, paste0("pars_", i, ".csv"))
      write.table(pars, file = fp, append = TRUE, sep = ",",
                  row.names = FALSE, col.names = FALSE)
      imps <- unlist(sapply(mdf@variables, FUN = function(y) {
        if(is(y, "irrelevant")) return(NA_real_)
        else if(y@all_obs) return(NA_real_)
        else return(y@imputations[s,,drop=TRUE])
      }))
      imps <- t(imps[!is.na(imps)])
      fp <- file.path(mdf@workpath, paste0("imps_", i, ".csv"))
      write.table(imps, file = fp, append = TRUE, sep = ",",
                  row.names = FALSE, col.names = FALSE)
    }
    Time.Elapsed <- proc.time() - ProcStart
    if(((Time.Elapsed)/60)[3] > max.minutes) {
      warning("'max.minutes' threshold exceeded")
      break
    }
  }
  if(((Time.Elapsed)/60)[3] > max.minutes) mdf@variables <- lapply(mdf@variables, .prune_missing_variable, s = s)
  
  if(verbose) message("Estimating models on completed data for chain ", i, "\n")

  mdf@variables <- lapply(mdf@variables, FUN = function(y) {
    if(!y@all_obs & !is(y, "irrelevant")) {
      model <- fit_model(y, mdf, s = s + 1, warn = TRUE)
      y@fitted <- fitted(model)
      if(!isS4(model)) model$x <- model$X <- model$y <- model$model <- NULL
      if(save_models) y@model <- model
    }
    else y@model <- NULL
    return(y)
  })
  
  mdf@done <- TRUE
  if(verbose) message("Done with chain ", i, "\n")
  return(mdf)
}

.mi_split <- function(i, y, data, verbose, s_start, s_end, ProcStart, max.minutes, parallel, save_models) {
  mdf <- y
  if(verbose) message("Chain ", i, "\n")
  data@priors <- mdf@priors
  for(s in s_start:s_end) {
    if(verbose) message("Chain ", i, " Iteration ", s, "\n")
    mdf <- fit_model(mdf, data, s = s, verbose = FALSE, warn = s == s_end)
    Time.Elapsed <- proc.time() - ProcStart
    if(((Time.Elapsed)/60)[3] > max.minutes) {
      warning("'max.minutes' threshold exceeded")
      break
    }
  }
  
  if(((Time.Elapsed)/60)[3] > max.minutes) mdf@variables <- lapply(mdf@variables, .prune_missing_variable, s = s)
  
  if(verbose) message("Estimating models on completed data for chain ", i, "\n")
  
  mdf@variables <- lapply(mdf@variables, FUN = function(y) {
    if(!y@all_obs & !is(y, "irrelevant")) {
      model <- fit_model(y, data, s = s + 1, warn = TRUE)
      y@fitted <- fitted(model)
      if(!isS4(model)) model$x <- model$X <- model$y <- model$model <- NULL
      if(save_models) y@model <- model
    }
    else y@model <- NULL
    return(y)
  })
  
  mdf@done <- TRUE
  if(verbose) message("Done with chain ", i, "\n")
  return(mdf)
}


setMethod("mi", signature(y = "missing_data.frame", model = "missing"), def =
  function(y, n.iter = 30,  n.chains = 4, max.minutes = Inf, seed = NA, verbose = TRUE, 
           save_models = FALSE, parallel = .Platform$OS.type != "Windows")
  {
    call <- match.call()
    if(!is.na(seed)) set.seed(seed)
    if(n.iter < 0) stop(message="number of iterations must be non-negative")
    ProcStart <- proc.time()
    s_start <- 0
    s_end <- n.iter

    Time.Elapsed <- proc.time() - ProcStart
    y@variables <- lapply(y@variables, FUN = function(x) {
      if(!x@all_obs & !is(x, "irrelevant")) {
        x@parameters  <- matrix(NA_real_, nrow = n.iter, ncol = 0)
        x@imputations <- matrix(NA_real_, nrow = n.iter, ncol = x@n_drawn)
        if(is(x, "semi-continuous")) {
          x@indicator@parameters  <- matrix(NA_real_, nrow = n.iter, ncol = 0)
          x@indicator@imputations <- matrix(NA_real_, nrow = n.iter, ncol = x@n_drawn)
        }
      }
      x@done <- TRUE
      return(x)
    })
    if(is(y, "allcategorical_missing_data.frame")) {
      y@latents@imputations <- matrix(NA_integer_, nrow = n.iter, ncol = nrow(y))
      y@latents@levels <- as.character(1:y@Hstar)
    }
    
    if(n.chains <= 0) return(y)
    
    if(is.logical(parallel) && parallel) {
      cores <- getOption("mc.cores", 2L)
      cl <- parallel::makeCluster(cores, outfile = "")
      on.exit(parallel::stopCluster(cl))
    }

    if(!parallel) {
      mdfs <- vector("list", n.chains)
      for(i in seq_along(mdfs)) {
        ProcStart <- proc.time()
        mdfs[[i]] <- .mi(i, y, verbose, s_start, s_end, ProcStart, 
                         max.minutes, parallel, save_models)
      }
    }
    else {
      mdfs <- parallel::parLapply(cl, X = as.list(1:n.chains),
                        fun = function(i) .mi(i, y, verbose, s_start, s_end,
                                              ProcStart, max.minutes, parallel, save_models))
    }
#     
#     else mdfs <- mclapply(as.list(1:n.chains),
#                           FUN = function(i) .mi(i, y, verbose, s_start, s_end,
#                                                 ProcStart, max.minutes, parallel, save_models))
    
    
    names(mdfs) <- paste("chain", 1:length(mdfs), sep = ":")
    
    object <- new("mi", 
                  call        = call,
                  data        = mdfs,
                  total_iters = as.integer(s_end))

    return(object)
  })

setMethod("mi", signature(y = "data.frame", model = "missing"), def =
  function(y, n.iter = 30,  n.chains = 4, max.minutes = Inf, seed = NA, verbose = TRUE, 
           save_models = FALSE, parallel = .Platform$OS.type != "Windows")
  {
    y <- as(y, "missing_data.frame")
    return(mi(y, n.iter = n.iter, n.chains = n.chains, max.minutes = max.minutes, 
              seed = seed, verbose = verbose, save_models = save_models, 
              parallel = parallel))
  })

setMethod("mi", signature(y = "matrix", model = "missing"), def =
  function(y, n.iter = 30,  n.chains = 4, max.minutes = Inf, seed = NA, verbose = TRUE,
           save_models = FALSE, parallel = .Platform$OS.type != "Windows")
  {
    y <- as(y, "missing_data.frame")
    return(mi(y, n.iter, n.chains, max.minutes, seed, verbose, save_models, parallel))
  })

setMethod("mi", signature(y = "mi", model = "missing"), 
          function(y, n.iter = 30,   max.minutes = Inf, seed = NA, verbose = TRUE, 
                   save_models = FALSE, parallel = .Platform$OS.type != "Windows")
              
          {
            call <- match.call()
            
            if(!is.na(seed)) set.seed(seed)
            
            if(n.iter < 1) stop(message="number of iterations must be at least 1")
            
            ProcStart <- proc.time()
            
            total_iters <- y@total_iters
            s_start <- sum(total_iters) + 1
            s_end <- s_start + n.iter - 1
            
            mdfs <- y@data
            n.chains <- length(mdfs)
            for(i in 1:n.chains) {
              y <- mdfs[[i]]
              if(TRUE) y@variables <- lapply(y@variables, FUN = function(x) {                
                if(x@all_obs & is(x, "irrelevant")) return(x)
                  x@imputations <- rbind(x@imputations, matrix(NA_integer_, n.iter, x@n_drawn))
                  x@parameters  <- rbind(x@parameters,  matrix(NA_real_, n.iter, ncol(x@parameters)))
                  if(is(x, "semi-continuous")) {
                    x@indicator@imputations <- rbind(x@indicator@imputations, matrix(NA_integer_, n.iter, x@indicator@n_drawn))
                    x@indicator@parameters  <- rbind(x@indicator@parameters, matrix(NA_real_, n.iter, ncol(x@indicator@parameters)))
                  }
                return(x)
              })
            }
            
            if(is.logical(parallel) && parallel) {
              cores <- getOption("mc.cores", 2L)
              cl <- parallel::makeCluster(cores, outfile = "")
              on.exit(parallel::stopCluster(cl))
            }

            if(!parallel) {
              mdfs <- vector("list", n.chains)
              for(i in seq_along(mdfs)) {
                ProcStart <- proc.time()
                mdfs[[i]] <- .mi(i, y, verbose, s_start, s_end, ProcStart, max.minutes, 
                                 parallel, save_models)
              }
            }
            else {
              mdfs <- parallel::parLapply(cl, as.list(1:n.chains),
                                fun = function(i) .mi(i, y, verbose, s_start, s_end,
                                                      ProcStart, max.minutes, parallel, save_models))
            }
#             else mdfs <- mclapply(as.list(1:n.chains),
#                                   FUN = function(i) .mi(i, y, verbose, s_start, s_end,
#                                                         ProcStart, max.minutes, parallel, save_models))
            
            object <- new("mi", 
                          call        = call,
                          data        = mdfs,
                          total_iters = as.integer(c(total_iters, n.iter)))
            return(object)
          })

setMethod("mi", signature(y = "missing_data.frame", model = "mi"), def =
  function(y, model, n.iter = sum(model@total_iters), max.minutes = 20, seed = NA, 
           verbose = TRUE, save_models = FALSE, 
           parallel = .Platform$OS.type != "Windows")
  {
    n.chains <- length(model)
    call <- match.call()
    if(!is.na(seed)) set.seed(seed)
    
    y <- mi(y, n.chains = 0L, n.iter = n.iter)
    
    ProcStart <- proc.time()
    s_start <- 0
    s_end <- n.iter

    if(is.logical(parallel) && parallel) {
      cores <- getOption("mc.cores", 2L)
      cl <- parallel::makeCluster(cores, outfile = "")
      on.exit(parallel::stopCluster(cl))
    }

    mdfs <- model@data
    if(!parallel) {
      for(i in seq_along(mdfs)) {
        ProcStart <- proc.time()
        mdfs[[i]] <- .mi_split(i, y, mdfs[[i]], verbose, s_start, s_end, ProcStart, 
                         max.minutes, parallel, save_models)
      }
    }
    else {
      mdfs <- parallel::parLapply(cl, as.list(1:n.chains),
                        fun = function(i) .mi_split(i, y, mdfs[[i]], verbose, s_start, s_end,
                                              ProcStart, max.minutes, parallel, save_models))
    }
#     else mdfs <- mclapply(as.list(1:n.chains),
#                           FUN = function(i) .mi_split(i, y, mdfs[[i]], verbose, s_start, s_end,
#                                                 ProcStart, max.minutes, parallel, save_models))
    
    names(mdfs) <- paste("chain", 1:length(mdfs), sep = ":")
    
    to_drop <- 1:ncol(model@data[[1]]@X)
    for(i in 1:n.chains) {
      model@data[[i]]@variables <- c(model@data[[i]]@variables, mdfs[[i]]@variables)
      model@data[[i]]@no_missing <- c(model@data[[i]]@no_missing, mdfs[[i]]@no_missing)
      # leave patterns as is I guess
      model@data[[i]]@DIM[2] <- model@data[[i]]@DIM[2] + mdfs[[i]]@DIM[2]
      model@data[[i]]@DIMNAMES[[2]] <- c(model@data[[i]]@DIMNAMES[[2]], mdfs[[i]]@DIMNAMES[[2]])
      mdfs[[i]]@index <- lapply(mdfs[[i]]@index, FUN = function(x) if(is.null(x)) x else to_drop)
      model@data[[i]]@index <- c(model@data[[i]]@index, mdfs[[i]]@index)
      model@data[[i]]@weights <- c(model@data[[i]]@weights, mdfs[[i]]@weights)
      model@data[[i]]@priors <- c(model@data[[i]]@priors, mdfs[[i]]@priors)
    }
    
    object <- new("mi", 
                  call        = call,
                  data        = model@data,
                  total_iters = as.integer(s_end))
    
    return(object)
    
  })

setMethod("mi", signature(y = "mdf_list", model = "missing"), def =
  function (y, ...) {
    out <- lapply(y, FUN = mi, ...)
    class(out) <- "mi_list"
    return(out)
  })

setMethod("mi", signature(y = "list", model = "missing"), def =
  function (y, ...) {
    if(!all(sapply(y, is, class2 = "mi"))) {
      stop("all elements of 'y' must be mi objects or missing_data.frame objects")
    }
    
    ## FIXME: should probably check that all the mi objects are based on the same missing_data.frame
    mdfs <- lapply(mi, FUN = function(x) return(x@data))
    
    object <- new("mi", 
                  call        = y[[1]]@call,
                  data        = mdfs,
                  total_iters = y[[1]]@total_iters)
    return(object)
  })

setMethod("mi", signature(y = "mdf_list", model = "missing"), 
          function (y, n.iter = 30,  n.chains = 4, max.minutes = Inf, seed = NA, verbose = TRUE, 
                    save_models = FALSE, parallel = .Platform$OS.type != "Windows")
          {
            
            out <- lapply(y, mi, n.iter = n.iter, n.chains = n.chains, max.minutes = max.minutes, 
                          seed = seed, verbose = verbose, save_models = save_models, parallel = parallel)
            class(out) <- "mi_list"
            return(out)
          })

setMethod("mi", signature(y = "mi_list", model = "missing"), def = 
  function (y, ...) {
    out <- lapply(y, FUN = mi, ...)
    class(out) <- "mi_list"
    return(out)
  })

setMethod("show", signature(object = "mi"), def =
  function(object) {
    cat("Object of class", class(object), "with", length(object@data), "chains, each with", 
        sum(object@total_iters), "iterations.\n")
    cat("Each chain is the evolution of an object of", class(object@data[[1]]), "class with",
        nrow(object@data[[1]]), "observations on", ncol(object@data[[1]]), "variables.\n")
    return(invisible(NULL))
  })

setMethod("show", signature(object = "mi_list"), def =
  function(object) {
    sapply(object, show)
    return(invisible(NULL))
  })

setMethod("summary", signature(object = "mi"), def =
  function(object) {
    mdf <- object@data[[1]]
    matrices <- complete(object, to_matrix = TRUE, include_missing = FALSE)
    chains <- length(matrices)
    matrices <- array(unlist(matrices), dim = c(dim(mdf), chains), dimnames = c(dimnames(mdf), NULL))
    out <- vector("list", ncol(mdf))
    names(out) <- colnames(mdf)
    for(i in seq_along(out)) {
      if(mdf@no_missing[i]) {
        if(is(mdf@variables[[i]], "categorical")) {
          mat <- table(matrices[,i,1])
          lev <- mdf@variables[[i]]@levels
          if(length(lev) && length(dim(mat)) > 1) colnames(mat) <- lev
        }
        else mat <- summary(matrices[,i,1])
        out[[i]] <- list(is_missing = "all values observed", observed = mat)
      }
      else if(is(mdf@variables[[i]], "categorical")) {
        mark <- is.na(mdf@variables[[i]])
        mat <- table(c(matrices[,i,]), rep(mark, times = chains))
        lev <- mdf@variables[[i]]@levels
        if(length(lev)) rownames(mat) <- lev
        colnames(mat) <- c("observed", "imputed")
        out[[i]] <- list(crosstab = mat) 
      }
      else {
        missing <- is.na(mdf@variables[[i]]@raw_data)
        out[[i]] <- list(is_missing = table(missing), imputed = summary(c(matrices[missing,i,])), 
                         observed = summary(c(matrices[!missing,i,])))
      }
    }
    return(out)
  })

setMethod("traceplot", signature(x = "mi"), def =
  function(x, ...) {
    traceplot(mi2BUGS, ...)
  })

setMethod("traceplot", signature(x = "mi_list"), def =
  function(x, ...) {
    traceplot(lapply(x, mi2BUGS, ...))
  })


## all the mi() methods below should return the missing_variable after imputing
## need to explicitly write out methods instead of doing poor man's S4
setMethod("mi", signature(y = "missing_variable", model = "ANY"), def = 
  function(y, model, ...) {
    stop("This method should not have been called. You need to define the relevant mi() S4 method")
  })

setMethod("mi", signature(y = "missing_variable", model = "missing"), def = 
  function(y) {
    if(y@n_drawn == 0) stop("'impute' should not have been called because there are no missing data")
    draws <- sample(y@data[y@which_obs], size = y@n_drawn, replace = TRUE)
    y@data[y@which_drawn] <- draws
    return(y)
  })

setMethod("mi", signature(y = "semi-continuous", model = "missing"), def = 
  function(y) {
    if(y@n_drawn == 0) stop("'impute' should not have been called because there are no missing data")
    y@indicator <- mi(y@indicator)
    draws <- sample(y@data[y@which_obs], size = y@n_drawn, replace = TRUE)
    if(is(y, "SC_proportion")) {
      n <- y@n_total
      if(is(y@indicator, "binary")) {
        mark <- which(complete(y@indicator, m = 0L)[y@which_miss] == 1)
        if(any(y@raw_data == 0, na.rm = TRUE)) draws[mark] <- .5 / n
        else                                   draws[mark] <- (n - .5) / n
      }
      else {
        mark <- which(complete(y@indicator, m = 0L)[y@which_miss] != 0)
        draws[mark] <- (draws[mark] * (n - 1) + .5) / n
      }
    }
    else if(is(y, "nonnegative-continuous")) {
      mark <- which(y@indicator@data[y@which_miss] == 1)
      if(length(mark)) draws[mark] <- y@transformation(rep(0, length(mark)))
    }
    else stop("FIXME: semi-continuous is not supported yet")
    y@data[y@which_drawn] <- draws
    return(y)
  })

# setMethod("mi", signature(y = "semi-continuous", model = "missing"), def = 
# function(y) {
#   if(y@n_drawn == 0) stop("'impute' should not have been called because there are no missing data")
# 
#   categories <- 1:(ncol(y@indicator@dummies) + 1)
#   draws <- sample(categories, size = y@n_drawn, replace = TRUE)
#   dummies <- t(sapply(draws, FUN = function(x) x == categories))[,-1,drop = FALSE]
#   y@indicator@dummies[y@which_drawn,] <- dummies
#   y@indicator@data[y@which_drawn] <- draws
# 
#   draws <- sample(y@data[y@which_obs], size = y@n_drawn, replace = TRUE)
#   if(is(y, "SC_proportion")) {
#     n <- y@n_total
#     if(is(y@indicator, "binary")) {
#       mark <- which(complete(y@indicator, m = 0L)[y@which_miss] == 1)
#       if(any(y@raw_data == 0, na.rm = TRUE)) draws[mark] <- .5 / n
#       else                                   draws[mark] <- (n - .5) / n
#     }
#     else {
#       mark <- which(complete(y@indicator, m = 0L)[y@which_miss] != 0)
#       draws[mark] <- (draws[mark] * (n - 1) + .5) / n
#     }
#   }
#   else if(is(y, "nonnegative-continuous")) {
#     mark <- which(y@indicator@data[y@which_miss] == 1)
#     if(length(mark)) draws[mark] <- y@transformation(rep(0, length(mark)))
#   }
# 
#   the_range <- range(y@data, na.rm = TRUE)
#   free <- y@data[y@which_obs]
#   free <- free[free != the_range[1] & free != the_range[2]]
#   draws <- sample(free, size = y@n_drawn, replace = TRUE)
#   y@data[y@which_drawn] <- draws
#   return(y)
# })

setMethod("mi", signature(y = "bounded-continuous", model = "missing"), def = 
  function(y) {
    if(y@n_drawn == 0) stop("'impute' should not have been called because there are no missing data")
    a <- if(length(y@lower) == 1) y@lower else y@lower[y@which_drawn]
    a <- ifelse(a == -Inf, min(y@data, na.rm = TRUE), a)
    a <- ifelse(a == Inf,  max(y@data, na.rm = TRUE), a)
    
    b <- if(length(y@upper) == 1) y@upper else y@upper[y@which_drawn]
    b <- ifelse(b == -Inf, min(y@data, na.rm = TRUE), b)
    b <- ifelse(b == Inf,  max(y@data, na.rm = TRUE), b)
    
    draws <- runif(y@n_drawn, min = a, max = b)
    y@data[y@which_drawn] <- draws
    return(y)
  })

setMethod("mi", signature(y = "categorical", model = "missing"), def = 
  function(y) {
    if(y@n_drawn == 0) stop("'impute' should not have been called because there are no missing data")
    draws <- sample(y@data[y@which_obs], size = y@n_drawn, replace = TRUE)
    y@data[y@which_drawn] <- draws
    return(y)
  })

.draw_parameters <-
  function(means, ev) {
    if(any(ev$values <= 0)) return(means)
    else return(means + (ev$vectors %*% (sqrt(ev$values) * rnorm(length(means))))[,1])
  }

.pmm <-
  function(y, eta, Sigma_inv = NULL, strata = NULL) {
    if(is(y, "unordered-categorical")) {
      if(is.null(Sigma_inv)) Sigma_inv <- .MPinverse(eta)
      MD <- mahalanobis(eta, colMeans(eta), Sigma_inv, inverted=TRUE)
      MD_observed <- MD[y@which_obs]
      y_observed <- y@data[y@which_obs]
      draws <- sapply(MD[y@which_drawn], FUN = function(x) {
        mark <- which.min(abs(MD_observed - x))
        drawmark <- c(y_observed[mark], mark)
        return(drawmark)
      })
    }
    else if(is(y, "grouped-binary")) {
      draws <- sapply(y@which_drawn, FUN = function(i) {
        which_same <- which(strata == strata[i])
        candidates <- intersect(which_same, y@which_obs)
        if(length(candidates) == 0) {
          msg <- paste(y@variable_name, ": must have some observed values in each group")
          stop(msg)
        }
        eta_can <- eta[candidates]
        y_can <- y@data[candidates]
        mark <- which.min(abs(eta_can - eta[i]))
        drawmark <- c(y_can[mark], mark)
        return(drawmark)
      })
    }
    else {
      eta_obs <-  eta[y@which_obs]
      y_obs <- y@data[y@which_obs]
      draws <- sapply(eta[y@which_drawn], FUN = function(x) {
        if(is.na(x)) return(NA_real_) # happens with semi-continuous
        mark <- which.min(abs(eta_obs - x))
        drawmark <- c(y_obs[mark], mark)
        return(drawmark)
      })
    }
    return(t(draws))
  }

setOldClass("polr")
setMethod("mi", signature(y = "ordered-categorical", model = "polr"), def = 
  function(y, model, s, ...) {
    if(y@n_drawn == 0) stop("'impute' should not have been called because there are no missing data")
    
    if(!is.element(y@imputation_method,  c("ppd", "pmm"))) badHessian <- FALSE
    else if(is.null(model$Hessian)) badHessian <- FALSE
    else if(!all(is.finite(model$Hessian))) badHessian <- TRUE
    else {
      means <- c(coef(model), model$zeta)
      ev <- eigen(vcov(model), symmetric = TRUE)
      badHessian <- any(ev$values <= 0)
      parameters <- .draw_parameters(means, ev)
      while(!badHessian && 
        any(diff(parameters[-(1:ncol(model$x))]) <= 0)) { # rejection sampling on cutpoints
        parameters <- .draw_parameters(means, ev)
      }
    }
    
    if(badHessian && y@imputation_method == "ppd") {
      warning(paste("predictive mean matching used for", y@variable_name, "on iteration", s,
                    "as a fallback due to Hessian error"))
      
      old_method <- y@imputation_method
      y@imputation_method <- "pmm"
      y <- mi(y, model, s, ...)
      y@imputation_method <- old_method
      return(y)
    }
    else if(y@imputation_method == "ppd") {
      eta <- as.vector(model$x[y@which_drawn,,drop=FALSE] %*% head(parameters, ncol(model$x)))
      pfun <- switch(y@family$link, logit = plogis, probit = pnorm, 
                     cloglog = function(q) exp(-exp(-q)), cauchit = pcauchy)
      zeta <- parameters[-(1:ncol(model$x))]
      draws <- sapply(eta, FUN = function(x) {
        which(rmultinom(1, 1, diff(c(0,pfun(zeta - x),1))) == 1)
      })
    }
    else if(y@imputation_method == "pmm") {
      parameters <- c(coef(model), model$zeta)
      eta <- model$x %*% parameters[1:ncol(model$x)]
      pmm <- .pmm(y, eta)
      draws <- pmm[,1]
      y@fitted[y@which_drawn,] <- y@fitted[y@which_obs,][pmm[,2],]
    }
    else if(y@imputation_method == "median") {
      predictions <- predict(model, type = "class")
      draws <- rep(floor(median(predictions[y@which_obs])), y@n_drawn)
    }
    else if(y@imputation_method == "mode") draws <- predict(model, type = "class")[y@which_drawn]
    else stop("'imputation_method' not recognized")
    
    y@data[y@which_drawn] <- draws
    y@imputations[s,] <- draws
    return(y)
  })

setOldClass("multinom")
setMethod("mi", signature(y = "unordered-categorical", model = "multinom"), def = 
  function(y, model, s, ...) {
    if(y@n_drawn == 0) stop("'impute' should not have been called because there are no missing data")
    ev <- eigen(vcov(model), symmetric = TRUE)
    parameters <- .draw_parameters(t(coef(model)), ev)
    if (ncol(model.matrix(model)) != nrow(parameters)) parameters <- t(parameters)
    eta <- model.matrix(model) %*% parameters
    if(y@imputation_method == "ppd") {
      exp_eta <- matrix(pmin(.Machine$double.xmax / ncol(eta), 
                             cbind(1, exp(eta[y@which_drawn,,drop = FALSE]))), 
                        ncol = ncol(eta) + 1)
      denom <- rowSums(exp_eta)
      Pr <- exp_eta / denom
      if (y@use_NA) {
        Pr <- Pr[,-1]/rowSums(Pr[,-1])
        badrows <- apply(is.na(Pr), 1, all)
        if(any(badrows)) {
          warning("Some rows of Pr are all 0 after dropping the missingness category")
          Pr[badrows,] <- 1/(ncol(Pr) - 1)
        }
      }
      draws <- apply(Pr, 1, FUN = function(p) which(rmultinom(1, 1, p) == 1))
    }
    else if(y@imputation_method == "pmm"){
    	pmm <- .pmm(y, eta)
    	draws <- pmm[,1]
    	y@fitted[y@which_drawn,,drop=FALSE] <- y@fitted[y@which_obs,,drop=FALSE][pmm[,2]]
    } 
    else if(y@imputation_method == "mode") draws <- predict(model, type = "class")[y@which_drawn]
    else stop("'imputation_method' not recognized")
    
    y@data[y@which_drawn] <- draws
    y@imputations[s,] <- draws
    return(y)
  })

setOldClass("RNL")
setMethod("mi", signature(y = "unordered-categorical", model = "RNL"), def = 
  function(y, model, s, ...) {
    if(y@n_drawn == 0) stop("'impute' should not have been called because there are no missing data")
    if(y@imputation_method == "ppd") { # imputating from the posterior predictive distribution
      Pr <- sapply(model, FUN = function(m) {
      	  ev <- eigen(vcov(m), symmetric = TRUE)
    	  parameters <- .draw_parameters(coef(m), ev)
    	  eta <- m$x[y@which_drawn,,drop=FALSE] %*% parameters
    	  pred <- m$family$linkinv(eta)
    	  return(pred)
  	  })
      if(y@use_NA) {
        Pr <- Pr[,-1]/rowSums(Pr[,-1])
        badrows <- apply(is.na(Pr), 1, all)
        if(any(badrows)) {
          warning("Some rows of Pr are all 0 after dropping the missingness category")
          Pr[badrows,] <- 1/(ncol(Pr) - 1)
        }
      }
      draws <- apply(Pr, 1, FUN = function(p) which(rmultinom(1, 1, p) == 1))
    }
    else if(y@imputation_method == "pmm") {
	  eta <- sapply(model, FUN = function(m) {
        ev <- eigen(vcov(m), symmetric = TRUE)
        parameters <- .draw_parameters(coef(m), ev)
        eta <- m$x %*% parameters
        return(eta)
      })
      pmm <- .pmm(y, eta)
      draws <- pmm[,1]
      y@fitted[y@which_drawn,,drop=FALSE] <- y@fitted[y@which_obs,,drop=FALSE][pmm[,2]]
    }
    else stop("only ppd and pmm are supported imputation methods in the RNL case")
    
    y@data[y@which_drawn] <- draws
    y@imputations[s,] <- draws
    return(y)
  })

setOldClass("glm")
setMethod("mi", signature(y = "binary", model = "glm"), def = 
  function(y, model, s, ...) {
    if(y@n_drawn == 0) stop("'impute' should not have been called because there are no missing data")
    if(y@imputation_method == "ppd") {
      ev <- eigen(vcov(model), symmetric = TRUE)
      parameters <- .draw_parameters(coef(model), ev)
      eta <- model$x[y@which_drawn,,drop=FALSE] %*% parameters
      pred <- model$family$linkinv(eta)
      draws <- rbinom(y@n_drawn, 1, pred) + 1L
    }
    else if(y@imputation_method == "pmm") {
      ev <- eigen(vcov(model), symmetric = TRUE)
      parameters <- .draw_parameters(coef(model), ev)
      eta <- model$x %*% parameters
      pmm <- .pmm(y, eta)
      draws <- pmm[,1]
      y@fitted[y@which_drawn] <- y@fitted[y@which_obs][pmm[,2]]
    }
    else if(y@imputation_method == "median") {
      predictions <- predict(model, type = "class")
      draws <- rep(floor(median(predictions[y@which_obs])), y@n_drawn)
    }
    else if(y@imputation_method == "mode") draws <- predict(model, type = "class")[y@which_drawn]
    else if(y@imputation_method == "mean") stop("'mean' is not a supported 'imputation_method' for binary variables")
    else if(y@imputation_method == "expectation") stop("'expectation' is not a supported 'imputation_method' for binary variables")
    else stop("'imputation_method' not recognized")
    
    draws <- as.integer(draws)
    y@data[y@which_drawn] <- draws
    y@imputations[s,] <- draws
    return(y)
  })

setOldClass("clogit")
setMethod("mi", signature(y = "grouped-binary", model = "clogit"), def = 
function(y, model, s, ...) {
  if(y@n_drawn == 0) stop("'impute' should not have been called because there are no missing data")
  # reconstruc the strata
  Terms <- model$terms
  temp <- untangle.specials(Terms, "strata")
  mf <- model.frame(model)
  strata <- strata(mf[, temp$vars], shortlabel = TRUE)
  if(y@imputation_method == "pmm") {
    ev <- eigen(vcov(model), symmetric = TRUE)
    parameters <- .draw_parameters(coef(model), ev)
    eta <- model$x %*% parameters
    draws <- .pmm(y, eta, strata = strata)[,1]
    #FIXME: haven't adjusted fitted values
  }
  else stop("only 'pmm' is supported for 'grouped-binary' variables")

  draws <- as.integer(draws)
  y@data[y@which_drawn] <- draws
  y@imputations[s,] <- draws
  return(y)
})

setMethod("mi", signature(y = "interval", model = "glm"), def = 
  function(y, model, s, ...) {
    stop("FIXME: write this method")
    if(y@n_drawn == 0) stop("'impute' should not have been called because there are no missing data")
    if(y@imputation_method == "ppd") {
      stop("FIXME")
    }
    else stop("only ppd is supported as an imputation method for interval variables")
    
    return(y)
  })

setMethod("mi", signature(y = "categorical", model = "matrix"), def = 
  function(y, model, s, ...) {
    if(y@n_drawn == 0) stop("'impute' should not have been called because there are no missing data")
    if(y@imputation_method != "ppd") stop("only ppd is supported in this case")
    if(nrow(model) != y@n_drawn) stop("matrix of probabilities has the wrong number of rows")
    draws <- apply(model, 1, FUN = function(p) which(rmultinom(1, 1, p) == 1))    
    y@data[y@which_drawn] <- draws
    y@imputations[s,] <- draws
    return(y)
  })

## helper function
.mi_continuous <-
  function(y, model) {
    if(y@n_drawn == 0) stop("'impute' should not have been called because there are no missing data")
    if(y@imputation_method == "ppd") {
      ev <- eigen(vcov(model), symmetric = TRUE)
      parameters <- .draw_parameters(coef(model), ev)
      if(model$family$family == "gaussian") {
        eta <- model$x[y@which_drawn,,drop=FALSE] %*% parameters
        pred <- model$family$linkinv(eta)
        if(is(y, "bounded-continuous")) {
          a <- if(length(y@lower) > 1) y@lower[y@which_drawn] else y@lower
          b <- if(length(y@upper) > 1) y@upper[y@which_drawn] else y@upper
          draws <- truncnorm::rtruncnorm(y@n_drawn, mean = pred, 
                                          sd = sqrt(model$dispersion), a = a, b = b)
        }
        else draws <- rnorm(y@n_drawn, pred, sqrt(model$dispersion))
      }
      else {
        eta <- model$x %*% parameters
        model$fitted <- model$family$linkinv(eta)
#         model$dispersion <- parameters@sigma^2
        draws <- y@family$sim(model, nsim = 1)[y@which_drawn]
      }
    }
    else if(y@imputation_method == "pmm") {
      ev <- eigen(vcov(model), symmetric = TRUE)
      parameters <- .draw_parameters(coef(model), ev)
      if(is(y, "semi-continuous")) {
        eta <- rep(NA_real_, y@n_total)
        mark <- complete(y@indicator, 0L) == 0
        eta[mark] <- model$x[mark,] %*% parameters
      }
      else eta <- model$x %*% parameters
      draws <- .pmm(y, eta)[,1]
      #FIXME: haven't adjusted fitted values using pmm for continuous
    }
    else if(y@imputation_method == "mean") {
      eta <- predict(model, type = "response")
      eta_observed <- eta[y@which_obs]
      eta_mean <- mean(eta_observed)
      draws <- rep(eta_mean, y@n_drawn)
    }
    else if(y@imputation_method == "median") {
      eta <- predict(model, type = "response")
      eta_observed <- eta[y@which_obs]
      eta_median <- median(eta_observed)
      draws <- rep(eta_median, y@n_drawn)
    }
    else if(y@imputation_method == "expectation") draws <- predict(model, type = "response")[y@which_drawn]
    else stop("'imputation_method' not recognized")
    return(draws)
  }

setMethod("mi", signature(y = "continuous", model = "glm"), def = 
  function(y, model, s, ...) {
    draws <- .mi_continuous(y, model)
    y@data[y@which_drawn] <- draws
    y@imputations[s,] <- draws
    return(y)
  })

# setMethod("mi", signature(y = "censored-continuous", model = "glm"), def = 
# function(y, model, s, ...) {
#   not_obs <- c(y@which_drawn, y@which_censored)
#   if(y@imputation_method == "ppd") {
#     parameters <- arm::sim(model, 1)
#     eta <- model$x[not_obs,,drop=FALSE] %*% parameters@coef[1,]
#     pred <- model$family$linkinv(eta)
#     draws <- rnorm(y@n_drawn, pred, parameters@sigma)
#   }
#   else if(y@imputation_method == "pmm") {
#     eta <- predict(model, type = "link")
#     eta_observed <- eta[y@which_obs]
#     y_observed <- y@data[y@which_obs]
#     draws <- sapply(eta[nob_obs], FUN = function(x) {
#                      mark <- which.min(abs(eta_observed - x))
#                      return(y_observed[mark])
#                    })
#   }
#   else if(y@imputation_method == "mean") {
#     eta <- predict(model, type = "response")
#     eta_observed <- eta[y@which_obs]
#     eta_mean <- mean(eta_observed)
#     draws <- rep(eta_mean, length(not_obs))
#   }
#   else if(y@imputation_method == "median") {
#     eta <- predict(model, type = "response")
#     eta_observed <- eta[y@which_obs]
#     eta_median <- median(eta_observed)
#     draws <- rep(floor(eta_median), length(not_obs))
#   }
#   else if(y@imputation_method == "expectation") draws <- predict(model, type = "response")[not_obs]
#   else stop("'imputation_method' not recognized")
# 
#   y@data[not_obs] <- draws
#   y@imputations[s,] <- draws
#   return(y)
# })

setMethod("mi", signature(y = "semi-continuous", model = "glm"), def = 
  function(y, model, s, ...) {
    stop("the semi-continuous mi() method should not have been called")
  })

setMethod("mi", signature(y = "nonnegative-continuous", model = "glm"), def = 
  function(y, model, s, ...) {
    draws <- .mi_continuous(y, model)
    # now account for the fact that some draws were determined to be 0 in step 1
    mark <- which(complete(y@indicator, 0L)[y@which_miss] == 1)
    if(length(mark)) draws[mark] <- y@transformation(rep(0, length(mark)))
    y@data[y@which_drawn] <- draws
    y@imputations[s,] <- draws
    return(y)
  })

## helper function
.mi_proportion <-
  function(y, model) {
    if(y@n_drawn == 0) stop("'impute' should not have been called because there are no missing data")
    
    if(!is.element(y@imputation_method,  c("ppd", "pmm"))) badHessian <- FALSE
    else if(is.null(model$vcov)) badHessian <- FALSE
    else if(!all(is.finite(model$vcov))) badHessian <- TRUE
    else {
      ev <- eigen(vcov(model), TRUE)
      badHessian <- any(ev$values <= 0)
      means <- coef(model)
      parameters <- .draw_parameters(means, ev)
#       while(!badHessian && parameters[length(parameters)] <= 0) {
#         parameters <- .draw_parameters(means, ev)
#       }
    }
    
    if(badHessian && y@imputation_method == "ppd") {
      warning(paste("predictive mean matching used for", y@variable_name,
                    "as a fallback due to Hessian error"))
      
      old_method <- y@imputation_method
      y@imputation_method <- "pmm"
      y <- mi(y, model)
      return(y@data[y@which_miss])
    }    
    else if(y@imputation_method == "ppd") {
      eta <- model$x[y@which_drawn,,drop=FALSE] %*% parameters[1:NCOL(model$x)]
      mu <- model$link$mean$linkinv(eta)
      phi <- model$link$precision$linkinv(parameters[length(parameters)]) ## FIXME: in the parameterized case
      shape1 <- mu * phi
      shape2 <- phi - shape1
      draws <- rbeta(y@n_drawn, shape1, shape2)
    }
    else if(y@imputation_method == "pmm") {
      eta <- model$x %*% parameters[-length(parameters)]
      draws <- .pmm(y, eta)[,1] #FIXME: haven't adjusted fitted values for pmm
    }
    else if(y@imputation_method == "mean") {
      mu <- predict(model)
      mu_observed <- mu[y@which_obs]
      mu_mean <- mean(mu_observed)
      draws <- rep(mu_mean, y@n_drawn)
    }
    else if(y@imputation_method == "median") {
      mu <- predict(model)
      mu_observed <- mu[y@which_obs]
      mu_median <- median(mu_observed)
      draws <- rep(mu_median, y@n_drawn)
    }
    else if(y@imputation_method == "expectation") draws <- predict(model)[y@which_drawn]
    else stop("'imputation_method' not recognized")
    return(draws)
  }

setOldClass("betareg")
setMethod("mi", signature(y = "proportion", model = "betareg"), def = 
  function(y, model, s, ...) {
    draws <- .mi_proportion(y, model)
    y@data[y@which_drawn] <- draws
    y@imputations[s,] <- draws
    return(y)
  })

setMethod("mi", signature(y = "proportion", model = "glm"), def = 
  function(y, model, s, ...) {
    draws <- .mi_continuous(y, model)
    y@data[y@which_drawn] <- draws
    y@imputations[s,] <- draws
    return(y)
  })

setMethod("mi", signature(y = "SC_proportion", model = "betareg"), def = 
  function(y, model, s, ...) {
    draws <- .mi_proportion(y, model)
    n <- y@n_total
    if(is(y@indicator, "binary")) {
      mark <- which(complete(y@indicator, 0L)[y@which_miss] == 1)
      if(any(y@raw_data == 0, na.rm = TRUE)) draws[mark] <- .5 / n
      else                                   draws[mark] <- (n - .5) / n
    }
    else {
      signs <- complete(y@indicator, 0L)[y@which_drawn]
      draws[signs < 0] <- .5 / n
      draws[signs > 0] <- (n - .5) / n
    }
    
    y@data[y@which_drawn] <- draws
    y@imputations[s,] <- draws
    return(y)
  })

## draw from overdispersed Poisson distribution
.rpois.od <- function(n, lambda, dispersion = 1) {
  if (dispersion <= 1) ans <- rpois(n, lambda)
  else {
    B <- 1/(dispersion-1)
    A <- lambda * B
    ans <- rnbinom(n, size= A , mu = lambda)
  }
  return(ans)
}

setMethod("mi", signature(y = "count", model = "glm"), def = 
  function(y, model, s, ...) {
    if(y@n_drawn == 0) stop("'impute' should not have been called because there are no missing data")
    if(y@imputation_method == "ppd") {
      ev <- eigen(vcov(model), symmetric = TRUE)
      parameters <- .draw_parameters(coef(model), ev)
      eta <- model$x[y@which_drawn,,drop=FALSE] %*% parameters
      pred <- model$family$linkinv(eta)
      draws <- .rpois.od(y@n_drawn, pred, model$dispersion)
    }
    else if(y@imputation_method == "pmm") {
      ev <- eigen(vcov(model), symmetric = TRUE)
      parameters <- .draw_parameters(coef(model), ev)
      eta <- model$x %*% parameters
      draws <- .pmm(y, eta)[,1] #FIXME: haven't adjusted fitted values for pmm
    }
    else if(y@imputation_method == "mean") {
      eta <- predict(model, type = "response")
      eta_observed <- eta[y@which_obs]
      eta_mean <- mean(eta_observed)
      draws <- rep(round(eta_mean), y@n_drawn)
    }
    else if(y@imputation_method == "median") {
      eta <- predict(model, type = "response")
      eta_observed <- eta[y@which_obs]
      eta_median <- median(eta_observed)
      draws <- rep(floor(eta_median), y@n_drawn)
    }
    else if(y@imputation_method == "expectation") draws <- round(predict(model, type = "response")[y@which_drawn])
    else stop("'imputation_method' not recognized")
    
    draws <- as.integer(draws)
    y@data[y@which_drawn] <- draws
    y@imputations[s,] <- draws
    return(y)
  })

setMethod("mi", signature(y = "irrelevant", model = "ANY"), def = 
  function(y, model, ...) {
    stop("The mi() method should not have been called on an 'irrelevant' variable")
  })

## FIXME: account for the other stuff at the bottom of the original mi.R file
