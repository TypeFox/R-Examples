## S4 StatModel object
RaschModel <- function(reltol = 1e-10, deriv = c("sum", "diff", "numeric"),
                       hessian = TRUE, maxit = 100L) {
  new("StatModel",
    capabilities = new("StatModelCapabilities"),
    name = "Rasch model",
    dpp = ModelEnvFormula,
    fit = function(object, weights = NULL, ...){

        ## extract response (there are no regressors)
        y <- object@get("response")

        ## call RaschModel.fit()
        z <- RaschModel.fit(y = y, weights = weights, reltol = reltol,
	  deriv = deriv, hessian = hessian, maxit = maxit)
        z$ModelEnv <- object
        z$addargs <- list(...)
        z
    }
  )
}

## methods needed for mob()
reweight.RaschModel <- function(object, weights, ...) {
     deriv <- if(is.null(object$deriv)) "sum" else object$deriv
     fit <- RaschModel(reltol = object$reltol, deriv = deriv)@fit
     do.call("fit", c(list(object = object$ModelEnv, weights = weights), object$addargs))
}

bread.RaschModel <- function(x, ...) x$vcov * x$n

estfun.RaschModel <- function(x, ...) {
  ## extract data and parameters of interest
  par <- x$coefficients
  esf <- x$elementary_symmetric_functions
  y <- x$data
  weights_orig <- weights(x)
  y <- y[weights_orig > 0, , drop = FALSE]
  weights <- weights_orig[weights_orig > 0]
  rs <- rowSums(y)
  
  ## analytical gradient
  if(!x$na) {
    agrad <- weights * (- y + esf[[2]][rs + 1, , drop = FALSE] / esf[[1]][rs + 1])[,-1, drop = FALSE]
  } else {
    ## set up return value
    n <- nrow(y)
    k <- ncol(y)
    agrad <- matrix(0, nrow = n, ncol = k)

    ## observed NA patterns
    na_patterns <- factor(apply(is.na(y), 1, function(z) paste(which(z), collapse = "\r")))

    ## loop over observed NA patterns	   
    for(i in seq_along(levels(na_patterns))) {
      ## parse NA pattern
      lev_i <- levels(na_patterns)[i]
      na_i <- which(na_patterns == lev_i)
      wi_i <- as.integer(strsplit(lev_i, "\r")[[1]])
      wi_i <- if(length(wi_i) < 1) 1:k else (1:k)[-wi_i]

      ## compute gradient per pattern
      esf_i <- esf[[i]]
      rs_i <- rowSums(y[na_i, wi_i, drop = FALSE])
      agrad[na_i, wi_i] <- weights[na_i] * (- y[na_i, wi_i, drop = FALSE] +
    	esf_i[[2]][rs_i + 1, , drop = FALSE] / esf_i[[1]][rs_i + 1])
    }

    agrad <- agrad[, -1, drop = FALSE]
  }

  ## collect and return
  grad <- matrix(0, ncol = length(par), nrow = length(weights_orig))
  grad[weights_orig > 0,] <- agrad
  return(grad)
}

