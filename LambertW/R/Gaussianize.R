#' @title Gaussianize matrix-like objects
#' 
#' @description
#' \code{Gaussianize} is probably the most useful function in this package. It
#'     works the same way as \code{\link[base]{scale}}, but instead of just
#'     centering and scaling the data, it actually \emph{Gaussianizes} the data
#'     (works well for unimodal data).  See Goerg (2011, 2016) and Examples.
#'
#' \strong{Important:} For multivariate input \code{X} it performs a column-wise
#'     Gaussianization (by simply calling \code{apply(X, 2, Gaussianize)}),
#'     which is only a marginal Gaussianization.  This does \emph{not} mean (and
#'     is in general definitely not the case) that the transformed data is then
#'     jointly Gaussian.
#' 
#' By default \code{Gaussianize} returns the \eqn{X \sim N(\mu_x, \sigma_x^2)}
#'     input, not the zero-mean, unit-variance \eqn{U \sim N(0, 1)} input.  Use
#'     \code{return.u = TRUE} to obtain \eqn{U}. 
#' @param data a numeric matrix-like object; either the data that should be
#'     Gaussianized; or the data that should ''DeGaussianized'' (\code{inverse =
#'     TRUE}), i.e., converted back to the original space.
#' 
#' @param type what type of non-normality: symmetric heavy-tails \code{"h"}
#'     (default), skewed heavy-tails \code{"hh"}, or just skewed \code{"s"}.
#' 
#' @param method what estimator should be used: \code{"MLE"} or \code{"IGMM"}.
#'     \code{"IGMM"} gives exactly Gaussian characteristics (kurtosis
#'     \eqn{\equiv} 3 for \code{"h"} or skewness \eqn{\equiv} 0 for \code{"s"}),
#'     \code{"MLE"} comes close to this. Default: \code{"IGMM"} since it is much
#'     faster than \code{"MLE"}.
#' 
#' @param return.tau.mat logical; if \code{TRUE} it also returns the estimated
#'     \eqn{\tau} parameters as a matrix (same number of columns as
#'     \code{data}).  This matrix can then be used to \code{Gaussianize} new
#'     data with pre-estimated \eqn{\tau}. It can also be used to
#'     ``DeGaussianize'' data by passing it as an argument (\code{tau.mat}) to
#'     \code{Gaussianize()} and set \code{inverse = TRUE}.
#' 
#' @param inverse logical; if \code{TRUE} it performs the inverse transformation
#'     using \code{tau.mat} to "DeGaussianize" the data back to the original
#'     space again.
#' 
#' @param tau.mat instead of estimating \eqn{\tau} from the data you can pass it
#'     as a matrix (usually obtained via \code{Gaussianize(..., return.tau.mat =
#'     TRUE)}). If \code{inverse = TRUE} it uses this \code{tau} matrix to
#'     ``DeGaussianize'' the data again.  This is useful to back-transform new
#'     data in the Gaussianized space, e.g., predictions or fits, back to the
#'     original space.
#' 
#' @param verbose logical; if \code{TRUE}, it prints out progress information in
#'     the console. Default: \code{FALSE}.
#' 
#' @param return.u logical; if \code{TRUE} it returns the zero-mean, unit
#'     variance Gaussian input.  If \code{FALSE} (default) it returns the input
#'     \eqn{X}.
#'
#' @param input.u optional; if you used \code{return.u = TRUE} in a previous
#'     step, and now you want to convert the data back to original space, then
#'     you have to pass it as \code{input.u}.  If you pass numeric data as
#'     \code{data}, \code{Gaussianize} assumes that \code{data} is the input
#'     corresponding to \eqn{X}, not \eqn{U}.
#'
#' @return 
#' numeric matrix-like object with same dimension/size as input \code{data}. 
#' If \code{inverse = FALSE} it is the Gaussianize matrix / vector; 
#' if \code{TRUE} it is the ``DeGaussianized'' matrix / vector.
#' 
#' The numeric parameters of mean, scale, and skewness/heavy-tail parameters
#'     that were used in the Gaussianizing transformation are returned as
#'     attributes of the output matrix: \code{'Gaussianized:mu'},
#'     \code{'Gaussianized:sigma'}, and for
#' 
#' \item{type = \code{"h"}:}{\code{'Gaussianized:delta'} & \code{'Gaussianized:alpha'},}
#' \item{type = \code{"hh"}:}{\code{'Gaussianized:delta_l'} and \code{'Gaussianized:delta_r'} & 
#' \code{'Gaussianized:alpha_l'} and \code{'Gaussianized:alpha_r'},}
#' \item{type = \code{"s"}:}{\code{'Gaussianized:gamma'}.}
#' 
#' They can also be returned as a separate matrix using \code{return.tau.mat =
#'     TRUE}. In this case \code{Gaussianize} returns a list with elements:
#'     \item{input}{Gaussianized input data \eqn{\boldsymbol x} (or
#'     \eqn{\boldsymbol u} if \code{return.u = TRUE}),} \item{tau.mat}{matrix
#'     with \eqn{\tau} estimates that we used to get \code{x}; has same number
#'     of columns as \code{x}, and 3, 5, or 6 rows (depending on
#'     \code{type='s'}, \code{'h'}, or \code{'hh'}).}
#'
#' @keywords univar multivariate
#' @export
#' @examples
#'
#' # Univariate example
#' set.seed(20)
#' y1 <- rcauchy(n = 100)
#' out <- Gaussianize(y1, return.tau.mat = TRUE)
#' x1 <- get_input(y1, c(out$tau.mat[, 1]))  # same as out$input
#' test_normality(out$input) # Gaussianized a Cauchy!
#' 
#' kStartFrom <- 20
#' y.cum.avg <- (cumsum(y1)/seq_along(y1))[-seq_len(kStartFrom)]
#' x.cum.avg <- (cumsum(x1)/seq_along(x1))[-seq_len(kStartFrom)]
#' 
#' plot(c((kStartFrom + 1): length(y1)), y.cum.avg, type="l" , lwd = 2, 
#'      main="CLT in practice", xlab = "n", 
#'      ylab="Cumulative sample average", 
#'      ylim = range(y.cum.avg, x.cum.avg))
#' lines(c((kStartFrom+1): length(y1)), x.cum.avg, col=2, lwd=2)
#' abline(h = 0)
#' grid()
#' legend("bottomright", c("Cauchy", "Gaussianize"), col = c(1, 2), 
#'        box.lty = 0, lwd = 2, lty = 1)
#' 
#' plot(x1, y1, xlab="Gaussian-like input", ylab = "Cauchy - output")
#' grid()
#'
#' # multivariate example
#' y2 <- 0.5 * y1 + rnorm(length(y1))
#' YY <- cbind(y1, y2)
#' plot(YY)
#' 
#' XX <- Gaussianize(YY, type = "hh")
#' plot(XX)
#' 
#' out <- Gaussianize(YY, type = "h", return.tau.mat = TRUE, 
#'                    verbose = TRUE, method = "IGMM")
#'                    
#' plot(out$input)
#' out$tau.mat
#' 
#' YY.hat <- Gaussianize(data = out$input, tau.mat = out$tau.mat,
#'                       inverse = TRUE)
#' plot(YY.hat[, 1], YY[, 1])
#' 

Gaussianize <- function(data = NULL, type = c("h", "hh", "s"), 
                        method = c("IGMM", "MLE"),
                        return.tau.mat = FALSE,
                        inverse = FALSE,
                        tau.mat = NULL,
                        verbose = FALSE,
                        return.u = FALSE,
                        input.u = NULL) {
  
  type <- match.arg(type)
  method <- match.arg(method)
  
  stopifnot(is.logical(return.tau.mat),
            is.logical(inverse),
            is.null(tau.mat) || is.numeric(tau.mat))
  
  data.is.u <- FALSE
  if (is.null(data)) {
    if (is.null(input.u)) {
      stop("You must either provide 'data' or 'input.u'.")
    } else {
      data <- input.u
      data.is.u <- TRUE
      overall.mean <- mean.default(input.u)
      if (lp_norm(overall.mean) > 0.1) {
        warning("It seems like your input data is not standardize correctly ",
                "(approx zero mean, variance 1). Pass it as 'data' argument or",
                "set 'return.u = TRUE' when you previously called ",
                "Gaussianize(data, return.u = TRUE, inverse = FALSE).")
      }
    }
  }

  if (is.null(ncol(data))) {
    data <- as.matrix(data)
  }
  if (is.null(colnames(data))) {
    colnames(data) <- paste0("Y", seq_len(ncol(data)))
  }
  # do the back-transformation; don't Gaussianize but do the inverse
  if (inverse) {
    if (verbose) {
      cat("Converting Gaussianized data back to original space.\n")
    }
    if (is.null(tau.mat)) {
      stop("If inverse = TRUE, tau.mat can not be NULL.")
    }
    if (!is.matrix(tau.mat)) {
      tau.mat <- as.matrix(tau.mat)
    }
    
    if (ncol(tau.mat) != ncol(data)) {
      stop("Number of columns in tau.mat must be the same as input data. \n",
           " Every column of tau.mat defines the transformation for ",
           "corresponding column of y (or input.u).")
    }
    if (data.is.u) {
      if (verbose) {
        cat("Input is assumed to zero-mean, unit variance U (not X). \n")
      }
      # multiply by sigma and add location
      data <- sweep(data, 2, tau.mat["sigma_x", ], "*")
      data <- sweep(data, 2, tau.mat["mu_x", ], "+")
    } 
    # initalize the output
    X <- data
    for (ii in seq_len(ncol(tau.mat))) {
      X[, ii] <- get_output(data[, ii], tau = unlist(tau.mat[, ii]))
    }
    # return the back-transformed data
    result <- X
  } else {  # do the Gaussianizing transformation
    if (verbose) {
      cat("Making data Gaussian.\n")
    }
    # if no tau.mat is there, then estimate it from data
    if (is.null(tau.mat)) {
      switch(method, 
             MLE = {
               .estimate_tau <- function(y) {
                 MLE_LambertW(y, type = type, distname = "normal")$tau
               }
             },
             IGMM = {
              .estimate_tau <- function(y) {
                IGMM(y, type = type)$tau
              }
            })
      if (verbose) {
        cat("Estimating parameters for", ncol(data), "variables ")
        start.time <- proc.time()
      }
      # estimate first column separately to get the names of the 'tau' vector
      tau.1 <- .estimate_tau(data[, 1])
      if (verbose) {
        elapsed.time <- proc.time() - start.time
        total.estimated.time <- elapsed.time[3] * ncol(data)
        names(total.estimated.time) <- "seconds"
        if (total.estimated.time > 60) {
          total.estimated.time <- c("minutes" = total.estimated.time / 60) 
        }
        cat("(will take ~", round(total.estimated.time, 1), 
            " ", names(total.estimated.time), "). \n", sep ="")
      } 
      
      if (ncol(data) > 1) {
        if (ncol(data) == 2) {
          tau.mat <- .estimate_tau(data[, -1])
        } else {
          tau.mat <- apply(data[, -1], 2, .estimate_tau)
        }
        tau.mat <- cbind(tau.1, tau.mat)
      } else {
        tau.mat <- cbind(tau.1)
      }
      rownames(tau.mat) <- names(tau.1)
      colnames(tau.mat) <- colnames(data)
      tau.mat <- round(tau.mat, 6)  # so that delta ~= 0 becomes actually 0
    } else { 
      if (verbose) {
        cat("Use provided tau.\n")
      }# use the provided tau
      stopifnot(is.matrix(tau.mat),
                ncol(tau.mat) == ncol(data))
      
      type.tmp <- tau2type(tau.mat[, 1])
      if (type.tmp != type) {
        warning("The tau.mat you provided suggests a type '", type.tmp,
                "' transformation.\n",
                "However, you specified 'type = ", type, "' as the argument.\n",
                "It will ignore this and use '", type.tmp,
                "' transformation instead.")
        type <- type.tmp
      }
    }
    # complete every tau
    tau.mat <- apply(tau.mat, 2, complete_tau, type)
    X <- data
    if (verbose) {
      cat("Converting to Gaussianized input data.\n")
    }
    for (ii in seq_len(ncol(tau.mat))) {
      X[, ii] <- get_input(data[, ii], tau = unlist(tau.mat[, ii]))
    }
    colnames(X) <- paste0(colnames(data), ".X")
    if (return.u) {
      if (verbose) {
        cat("Return zero-mean, unit variance U (not X). \n")
      }
      X <- sweep(X, 2, tau.mat["mu_x", ], "-")
      X <- sweep(X, 2, tau.mat["sigma_x", ], "/")
      colnames(X) <- paste0(colnames(data), ".U")
    }
    
    attr(X, "Gaussianized:mu") <- tau.mat["mu_x", ]
    attr(X, "Gaussianized:sigma") <- tau.mat["sigma_x", ]
    
    switch(type,
           s = {
             attr(X, "Gaussianized:gamma") <- tau.mat["gamma", ]
           },
           h = {
             attr(X, "Gaussianized:delta") <- tau.mat["delta", ]
             attr(X, "Gaussianized:alpha") <- tau.mat["alpha", ]
           },
           hh = {
             attr(X, "Gaussianized:delta_l") <- tau.mat["delta_l", ]
             attr(X, "Gaussianized:delta_r") <- tau.mat["delta_r", ]
             attr(X, "Gaussianized:alpha_l") <- tau.mat["alpha_l", ]
             attr(X, "Gaussianized:alpha_r") <- tau.mat["alpha_r", ]
           }
    )
    
    
    if (return.tau.mat) {
      result <- list(input = X,
                     tau.mat = tau.mat)
    } else {
      result <- X
    }
  }
  if (verbose) {
    cat("Done! \n\n")
  }
  return(result)
}
