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


pool <-
  function(formula, data, m = NULL, FUN = NULL, ...) {
    if(is.list(data)) {
      if(all(sapply(data, is, class2 = "mi"))) {
        if(is.null(m)) m <- length(data[[1]]@data)
        else m <- as.integer(m)
        dfs <- complete(data[[1]], m = m, to_matrix = FALSE)
        l <- length(data)
        if(l > 1) for(i in 2:l) {
          temp <- complete(data[[l]], m = m, to_matrix = FALSE)
          for(j in seq_along(temp)) dfs[[j]] <- rbind(dfs[[j]], temp[[j]])
        }
        data <- data[[1]]
      }
      else if(all(sapply(data, is.data.frame))) {
        dfs <- data
        m <- length(dfs) 
      }
      else {
        stop("if 'data' is a list it must be a list of mi objects or data.frames")
      }
    }
    else if(is(data, "mi")) {
      if(is.null(m)) m <- length(data@data)
      else m <- as.integer(m)
      dfs <- complete(data, m = m, to_matrix = FALSE)
    }
    if(!is(formula, "formula")) stop("'formula' must be a formula")
    
    dots <- list(...)
    
    if(is.null(FUN)) {
      if(!is(data, "mi")) stop("if 'data' is not of class 'mi', 'FUN' must be specified")
      yname <- as.character(formula)[2]
      if(!(yname %in% colnames(data@data[[1]]))) {
        stop(paste("no variable called", yname, "possibly due to typo or transformation,",
                   "in which case you need to specify 'FUN' explicitly"))
      }
      else y <- data@data[[1]]@variables[[yname]]
      
      if(!is.method_in_mi("fit_model", y = class(y), data = class(data@data[[1]]))) {
        stop(paste(yname, "seems to have a user-defined 'fit_model' method,",
                   "in which case 'FUN' must be specified explicitly"))
      }
      
      if(is(y, "unordered-categorical")) {
        FUN <- nnet::multinom
        fit <- "multinom"
      }
      else if(is(y, "binary") | is(y, "count") | is(y, "continuous")) {
        FUN <- arm::bayesglm
        fit <- "bayesglm"
        if(!("family" %in% names(dots))) dots$family <- y@family
      }
      else if(is(y, "interval")) {
        FUN <- survival::survreg
        fit <- "survreg"
      }
      else if(is(y, "ordered-categorical")) {
        FUN <- arm::bayespolr
        fit <- "bayespolr"
        if(!("method" %in% names(dots))) dots$method <- if(y@family$link == "logit") "logistic" else y@family$link
      }
    }
    else if(!is(FUN, "function")) stop("'FUN' must be a function or NULL")
    else fit <- deparse(substitute(FUN))
    
    models <- lapply(dfs, FUN = function(d) {
      dots$data <- d
      dots$formula <- formula
      do.call(FUN, args = dots)
    })
    
    summaries <- lapply(models, summary)
    pooled_summary <- summaries[[1]]
    if(is.list(pooled_summary)) for(i in seq_along(pooled_summary)) {
      if(is.numeric(pooled_summary[[i]])) {
        num <- lapply(summaries, FUN = function(x) x[[i]])
        if(is.matrix(pooled_summary[[i]])) {
          mat <- pooled_summary[[i]]
          arr <- array(unlist(num), dim = c(dim(mat), m))
          arr <- apply(arr, 1:2, mean)
          colnames(arr) <- colnames(mat)
          rownames(arr) <- rownames(mat)
          pooled_summary[[i]] <- arr
        }
        else if(length(pooled_summary[[i]]) > 1) {
          arr <- rowMeans(matrix(unlist(num), ncol = m))
          names(arr) <- names(pooled_summary[[i]])
          pooled_summary[[i]] <- arr
        }
        else pooled_summary[[i]] <- mean(unlist(num))
      }
    }
    else {
      pooled_summary <- list()
      warning("could not construct pooled_summary")
    }
    
    coefs <- sapply(models, get_parameters)
    variances <- sapply(models, FUN = function(x) diag(vcov(x)))
    W <- rowMeans(variances)
    B <- apply(coefs, 1, var) 
    ses <- sqrt(W  + B * (1 + 1/m))
    
    if(is(pooled_summary, "summary.glm") | is(pooled_summary, "summary.polr")) {
      pooled_summary$call <- match.call()
      pooled_summary$coefficients[,1:2] <- cbind(rowMeans(coefs), ses)
    }
    else if(is(pooled_summary, "summary.multinom")) {
      pooled_summary$call <- match.call()
      pooled_summary$coefficients <- cbind(coef = rowMeans(coefs), ses, z = NA_real_, p = NA_real_)
    }
    else warning("pooled_summary is probably bogus")
    
    if(ncol(pooled_summary$coefficients) >= 3) {
      if(colnames(pooled_summary$coefficients)[3] == "t value") {
        pooled_summary$coefficients[,3] <- tvalue <- pooled_summary$coefficients[,1] / ses 
        if(TRUE) {
          gamma <- (1 + 1/m) * B / ses^2
          df.r <- pooled_summary$df.residual
          v <- (m - 1) * (1 + m/(m + 1) * W / B)^2
          v_obs <- (1 - gamma) * (df.r + 1) / (df.r + 3) * df.r
          df.star <- 1/(1/v + 1/v_obs)
          if(ncol(pooled_summary$coefficients) == 4) {
            pooled_summary$coefficients[,4] <- 2 * pt(-abs(tvalue), df.star)
          }
          else pooled_summary$coefficients <- cbind(pooled_summary$coefficients, 
                                                    "p-value" = 2 * pt(-abs(tvalue), df.star))
        }
      }
      else {
        pooled_summary$coefficients[,3] <- zvalue <- pooled_summary$coefficients[,1] / ses
        if(ncol(pooled_summary$coefficients) == 4) {
          pooled_summary$coefficients[,4] <- 2 * pnorm(-abs(zvalue))
        }
        else pooled_summary$coefficients <- cbind(pooled_summary$coefficients,
                                                  "p-value" = 2 * pnorm(-abs(zvalue)) )
      }
    }
    kall <- match.call()
    kall[1] <- call(fit)
    out <- new("pooled", formula = formula, fit = fit, models = models,
               coefficients = rowMeans(coefs), ses = ses, pooled_summary = pooled_summary, call = kall)
    return(out)
  }

setMethod("display", signature(object = "pooled"), def =
  function(object, digits = 2, ...) {
    call <- object@call
    summ <- summary(object)
    coef <- object@pooled_summary$coefficients[,1:2]
    colnames(coef) <- c("coef.est", "coef.se")
    n <- summ$df.residual
    k <- summ$df[1]
    k.intercepts <- length(summ$zeta)
    
    print(call)
    pfround(coef, digits)
    if(k.intercepts > 0) {
      cat(paste("n = ", n, ", k = ", k, " (including ", k.intercepts, 
                " intercepts)\nresidual deviance = ", fround(summ$deviance, 1),
                ", null deviance is not computed by polr", "\n", sep = ""))
      return(invisible(NULL))
    }
    cat(paste("n = ", n, ", k = ", k, "\nresidual deviance = ", 
              fround(summ$deviance, 1), ", null deviance = ", fround(summ$null.deviance, 1), 
              " (difference = ", fround(summ$null.deviance - summ$deviance, 1), ")", "\n", sep = ""))
    dispersion <- summ$dispersion 
    if (dispersion != 1) {
      cat(paste("overdispersion parameter = ", fround(dispersion, 1), "\n", sep = ""))
      if (summ$family$family == "gaussian") {
        cat(paste("residual sd is sqrt(overdispersion) = ", 
                  fround(sqrt(dispersion), digits), "\n", sep = ""))
      }
    }
    return(invisible(NULL))
  })

setMethod("show", signature(object = "pooled"), def =
  function(object) {
    display(object)
    return(invisible(NULL))
  })

setMethod("summary", signature(object = "pooled"), def = 
  function(object, ...) {
    return(object@pooled_summary)
  })

setMethod("coef", signature(object = "pooled"), def = 
  function(object, ...) {
    return(object@coefficients)
  })

setMethod("vcov", signature(object = "pooled"), def = 
  function(object, ...) {
    return(object@vcov)
  })

setMethod("residuals", signature(object = "pooled"), def = 
  function(object, ...) {
    return(rowMeans(sapply(object@models, residuals)))
  })

setMethod("fitted", signature(object = "pooled"), def = 
  function(object, ...) {
    return(rowMeans(sapply(object@models, fitted, ...)))
  })

