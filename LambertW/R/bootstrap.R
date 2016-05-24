#' @title Bootstrap Lambert W x F estimates
#' 
#' @description
#' Analyzes the Lambert W x F for a given dataset based on bootstrapping. Depends
#' on the \pkg{boot} package and returns a \code{"boot"} object.
#' 
#' @export
#' @param object an object of class \code{"LambertW_fit"}; usually
#' output of \code{\link{IGMM}} or \code{\link{MLE_LambertW}}.
#' @param sample.size sample size of the bootstrap.  By default, equal to the 
#' original data length.
#' @param R number of replicates for the bootstrap. See 
#' \code{\link[boot]{boot}} for details.
#' @param ... additional arguments passed to \code{\link{boot}}.
#' @return
#' An object of class \code{"boot"} representing the bootstrap 
#' analysis of \eqn{\hat{\theta}} (or \eqn{\hat{\tau}}) of 
#' an Lambert W x F estimator (\code{LambertW_fit}).
#' 
#' @examples
#' \dontrun{
#' yy <- rLambertW(n = 1000, theta = list(delta = c(0.1), beta = c(2, 1)), 
#'                 distname = "normal")
#' mod.igmm <- IGMM(yy, type = "h")
#' boot.est <- bootstrap(mod.igmm, R = 100) 
#' # use summary and plot from 'boot' pkg
#' plot(boot.est, 3)
#' summary(boot.est)
#' }


bootstrap <- function(object, ...) {
  UseMethod("bootstrap")
}

#' @rdname bootstrap
#' @export
#' 
bootstrap.LambertW_fit <- function(object, 
                                   sample.size = length(object$data),
                                   R = 100,
                                   ...) {

  stopifnot(inherits(object, "LambertW_fit"),
            is.numeric(sample.size),
            length(sample.size) == 1,
            sample.size > 0)
  
  if (!requireNamespace("boot", quietly = TRUE)) {
    stop("The 'boot' package is not available.  Please install it",
         "to use this function.") 
  }
  
  num.obs <- length(object$data)
  if (sample.size > num.obs) {
    stop("Sample size for bootstrapping must be less or equal to ",
         "original data.\n Currently, bootstrap sample size was set ",
         "to ", sample.size, ", but data has only ", num.obs, " observations.")
  }
    
  all.args <- as.list(object$call[-1]) 
  fct.name <- as.list(object$call)[[1]]
  fct.name <- as.character(object$call)[1]
  new.data.arg <- "new.data"

  if (grepl("::", fct.name)) {
    warning("Stripped of package name prefix from function name. ",
            "If this does not work, please call '",
            fct.name, "' without the 'LambertW::' prefix and call ",
            "bootstrap(", deparse(substitute(object)), ") again.")
    fct.name <- tail(strsplit(fct.name, ":")[[1]], 1)
  }
  .auxEst <- function(data, indices = seq_len(length(data))) {
    
    # take the first sample.size samples from the index; this only
    # works for sample.size <= length(indices)
    indices <- head(indices, sample.size)
    
    pe.boot <- new.env()
    assign(new.data.arg, data[indices], envir = pe.boot)
    # relabel to 'new.data.arg'
    all.args[["y"]] <- as.symbol(new.data.arg)
    
    mod.est <- do.call(as.character(fct.name), all.args,
                       envir = pe.boot)
    if (object$method == "IGMM") {
      est.hat <- mod.est$tau
    } else if (object$method == "MLE") {
      est.hat <- mod.est$theta
      for (nn in names(mod.est$theta.fixed)) {
        est.hat[[nn]] <- NULL
      }
      est.hat <- flatten_theta(est.hat)
    } else {
      stop("Bootstrap for method '", object$method, 
           "' is not implemented yet.")
    }
    names(est.hat) <- gsub("beta.", "", names(est.hat))
    return(est.hat)
  }
  out <- boot::boot(data = object$data,
                    statistic = .auxEst,
                    R = R, ...)
  return(out)
}
