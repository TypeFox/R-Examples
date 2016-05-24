#' @importFrom stats lm.fit qnorm rnorm
#' @importFrom utils tail
#' @import Rcpp
#' @useDynLib SLOPE
NULL

#' SLOPE: Sorted L1 Penalized Estimation
#'
#' Performs variable selection using SLOPE (Sorted L1 Penalized Estimation).
#' Given a design matrix \eqn{X} and a response vector \eqn{y}, find the
#' coefficient vector \eqn{\beta} minimizing
#' \deqn{\frac{1}{2} \Vert X\beta - y \Vert_2^2 +
#'       \sigma \cdot \sum_{i=1}^p \lambda_i |\beta|_{(i)},}
#' where the \eqn{\lambda} sequence is chosen to control the false discovery
#' rate associated with nonzero components of \eqn{\beta}.
#'
#' @param X the \eqn{n}-by-\eqn{p} design matrix
#' @param y response vector of length \eqn{n}
#' @param fdr target FDR (false discovery rate)
#' @param lambda specifcation of \eqn{\lambda}, either one of "bhq" or "gaussian",
#'  or a vector of length \eqn{p}, sorted in decreasing order
#'  (see \code{\link{create_lambda}})
#' @param sigma noise level. If omitted, estimated from the data (see Details).
#' @param normalize whether to center the input data and re-scale the columns
#'  of the design matrix to have unit norm. Do not disable this unless you are
#'  certain that your data is appropriately pre-processed.
#' @param solver which SLOPE solver to use (see Details)
#' @param ... additional arguments to pass to the solver (see the relevant
#'  solver)
#'
#' @return An object of class \code{SLOPE.result}. This object is a list
#'  containing at least the following components:
#'  \item{lambda}{the \eqn{\lambda} sequence used}
#'  \item{lambda_method}{method of \eqn{\lambda} construction
#'   ("bhq", "gaussian", or "user")}
#'  \item{sigma}{(sequence of) noise level(s) used}
#'  \item{beta}{optimized coefficient vector \eqn{\beta}}
#'  \item{selected}{selected variables
#'    (variables \eqn{i} with \eqn{\beta_i > 0})}
#'
#' @details At present, two solvers for the SLOPE problem are supported. By
#' default, we use \code{\link{SLOPE_solver}}, which is mostly written in R but
#' uses a fast prox implemented in C. If you have MATLAB installed, it is also
#' possible to use the TFOCS solver for SLOPE. This requires the MATLAB package
#' TFOCS and the R package \code{R.matlab}.
#'
#' If the noise level is unknown, it is estimated from the data using one of
#' two methods. When \eqn{n} is large enough compared to \eqn{p}, the classical
#' unbiased estimate of \eqn{\sigma^2} is used. Otherwise, the
#' \emph{iterative SLOPE} algorithm is used, in which a decreasing sequence of
#' \eqn{\sigma^2} estimates is used until the set of selected variables
#' stabilizes. For details, see Section 3.2.3 of the SLOPE paper.
#'
#' @seealso \code{\link{SLOPE_solver}}
#' @export
SLOPE <- function(X, y, fdr = 0.20, lambda = 'gaussian', sigma = NULL,
                  normalize = TRUE, solver = c('default','matlab'), ...) {
  # Validate input types.
  if (is.data.frame(X)) {
    X.names = names(X)
    X = as.matrix(X, rownames.force = F)
  } else if (is.matrix(X))
    X.names = colnames(X)
  else
    stop('Input X must be a matrix or data frame')
  n = nrow(X); p = ncol(X)
  y = as.numeric(y)

  # Normalize input, if necessary.
  if (normalize) {
    X = normalize(X)
    y = y - mean(y)
  }

  # Create lambda sequence.
  if (is.character(lambda)) {
    lambda_method = lambda
    lambda = create_lambda(n, p, fdr, lambda)
  } else {
    lambda_method = 'user'
    lambda = as.numeric(lambda)
  }

  # Validate input constraints.
  stopifnot(length(y) == n, length(lambda) == p)
  if (is.unsorted(rev(lambda)))
    stop('Lambda sequence must be non-increasing');

  # Estimate the noise level, if possible.
  if (is.null(sigma) && n >= p + 30)
     sigma = estimate_noise(X, y)

  # Run the solver, iteratively if necessary.
  if (is.null(sigma)) {
    # Run Algorithm 5 of Section 3.2.3.
    selected = NULL
    repeat {
      selected.prev = selected
      sigma = c(sigma, estimate_noise(X[,selected,drop=F], y))
      result = SLOPE_solver_call(solver, X, y, tail(sigma,1) * lambda)
      selected = result$selected
      if (identical(selected, selected.prev))
        break
      if (length(selected)+1 >= n)
        stop('Selected >= n-1 variables. Cannot estimate variance.')
    }
  } else {
    result = SLOPE_solver_call(solver, X, y, sigma * lambda)
  }

  # Package up the results.
  beta = result$beta
  selected = result$selected
  if (!is.null(X.names))
    names(selected) = X.names[selected]

  structure(list(call = match.call(),
                 lambda = lambda,
                 lambda_method = lambda_method,
                 sigma = sigma,
                 beta = beta,
                 selected = selected),
            class = 'SLOPE.result')
}

# Helper function for invoking the SLOPE solver.
SLOPE_solver_call <- function(solver = c('default','matlab'), ...) {
  solver_name = match.arg(solver)
  if (solver_name == 'default') {
    result = SLOPE_solver(...)
    beta = result$x
    tol = 0 # Our solver sets un-selected beta's to *exactly* zero.
  } else if (solver_name == 'matlab') {
    beta = SLOPE_solver_matlab(...)
    tol = 1e-5 # FIXME: How to choose this?
  }
  selected = which(abs(beta) > tol)
  list(beta = beta, selected = selected)
}

# Print method for SLOPE results.
#' @export
#' @keywords internal
print.SLOPE.result <- function(x, ...) {
  result = x
  cat("\nCall:\n", paste(deparse(result$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  selected = result$selected
  if (length(selected)) {
    beta.selected <- result$beta[selected]
    if (is.null(names(selected)))
      names(beta.selected) <- as.character(selected)
    else
      names(beta.selected) <- names(selected)
    cat("Selected", length(beta.selected), "variables with coefficients:\n")
    print.default(beta.selected, print.gap = 2L, quote = FALSE)
  } else {
    cat("No selected variables\n")
  }
  cat("\n")
  invisible(result)
}