#' Sorted L1 solver
#'
#' Solves the sorted L1 penalized regression problem: given a matrix \eqn{A},
#' a vector \eqn{b}, and a decreasing vector \eqn{\lambda}, find the vector
#' \eqn{x} minimizing
#' \deqn{\frac{1}{2}\Vert Ax - b \Vert_2^2 +
#'       \sum_{i=1}^p \lambda_i |x|_{(i)}.}
#'
#' @param A an \eqn{n}-by-\eqn{p} matrix
#' @param b vector of length \eqn{n}
#' @param lambda vector of length \eqn{p}, sorted in decreasing order
#' @param initial initial guess for \eqn{x}
#' @param prox function that computes the sorted L1 prox
#' @param max_iter maximum number of iterations in the gradient descent
#' @param grad_iter number of iterations between gradient updates
#' @param opt_iter number of iterations between checks for optimality
#' @param tol_infeas tolerance for infeasibility
#' @param tol_rel_gap tolerance for relative gap between primal and dual
#'  problems
#'
#' @return An object of class \code{SLOPE_solver.result}. This object is a list
#'  containing at least the following components:
#'  \item{x}{solution vector \eqn{x}}
#'  \item{optimal}{logical: whether the solution is optimal}
#'  \item{iter}{number of iterations}
#'
#' @details This optimization problem is convex and is solved using an
#' accelerated proximal gradient descent method.
#'
#' @export
# Adapted from Adlas.R in original SLOPE distribution.
SLOPE_solver <- function(A, b, lambda, initial = NULL, prox = prox_sorted_L1,
                         max_iter = 10000, grad_iter = 20, opt_iter = 1,
                         tol_infeas = 1e-6, tol_rel_gap = 1e-6) {
  # Get problem dimension.
  n = ncol(A)

  # Get initial lower bound on the Lipschitz constant.
  x = with_seed(0, rnorm(n))
  x = x / sqrt(sum(x^2))
  x = t(A) %*% (A %*% x)
  L = sqrt(sum(x^2))

  # Initialize parameters and iterates.
  x.init = if (is.null(initial)) rep(0,n) else initial
  t      = 1
  eta    = 2
  x      = x.init
  y      = x
  Ax     = A %*% x
  f.prev = Inf
  iter   = 0
  optimal = FALSE

  # Main loop.
  repeat {
    # Compute the gradient at f(y).
    if ((iter %% grad_iter) == 0) # Includes first iterations
      r = (A %*% y) - b
    else
      r = (Ax + ((t.prev - 1) / t) * (Ax - Ax.prev)) - b
    g = t(A) %*% r
    f = as.double(crossprod(r)) / 2

    # Increment iteration count.
    iter = iter + 1

    # Check optimality conditions.
    if ((iter %% opt_iter) == 0) {
      # Compute 'dual', check infeasibility and gap.
      gs     = sort(abs(g), decreasing=TRUE)
      ys     = sort(abs(y), decreasing=TRUE)
      infeas = max(max(cumsum(gs-lambda)),0)

      # Compute primal and dual objective.
      obj_primal =  f + as.double(crossprod(lambda,ys))
      obj_dual   = -f - as.double(crossprod(r,b))

      # Check primal-dual gap.
      if ((abs(obj_primal - obj_dual)/max(1,obj_primal) < tol_rel_gap) &&
            (infeas < tol_infeas * lambda[[1]]))
        optimal = TRUE;
    }

    # Stopping criteria.
    if (optimal || (iter >= max_iter))
      break;

    # Store copies of previous values.
    Ax.prev = Ax
    x.prev = x; f.prev = f; t.prev = t

    # Lipschitz search.
    repeat {
      # Compute prox mapping.
      x = prox(y - (1/L)*g, lambda/L)
      d = x - y

      Ax = A %*% x
      r  = Ax-b
      f  = as.double(crossprod(r))/2
      q  = f.prev + as.double(crossprod(d,g)) + (L/2)*as.double(crossprod(d))

      if (q < f*(1-1e-12))
        L = L * eta
      else
        break
    }

    # Update.
    t <- (1 + sqrt(1 + 4*t^2)) / 2
    y <- x + ((t.prev - 1) / t) * (x - x.prev)
  }
  if (!optimal)
    warning('SLOPE solver reached iteration limit')

  # Package up the results.
  structure(list(x = as.vector(y),
                 optimal = optimal,
                 iter = iter,
                 infeas = infeas,
                 obj_primal = obj_primal,
                 obj_dual = obj_dual,
                 lipschitz = L),
            class = 'SLOPE_solver.result')
}

#' Interface to MATLAB sorted L1 solver
#'
#' Uses the \code{R.matlab} package to invoke the TFOCS sorted L1 solver.
#' See \code{\link{SLOPE_solver}} for a description of the sorted L1
#' optimization problem.
#'
#' @param A an \eqn{n}-by-\eqn{p} matrix
#' @param b vector of length \eqn{n}
#' @param lambda vector of length \eqn{p}, sorted in decreasing order
#' @param initial initial guess for \eqn{x}
#' @param matlab MATLAB client object (instance of class \code{Matlab}).
#'  If supplied, the client object should be connected to a running MATLAB
#'  server. If not supplied, a MATLAB server will be started (requires 'matlab'
#'  to be on the PATH).
#'
#' @return The solution vector \eqn{x}.
#'
#' @seealso \code{\link{SLOPE_solver}}
#' @export
#' @keywords internal
SLOPE_solver_matlab <- function(A, b, lambda, initial=numeric(), matlab=NULL) {
  loadNamespace('R.matlab')

  if (is.null(matlab)) {
    R.matlab::Matlab$startServer()
    matlab = R.matlab::Matlab()
    open(matlab)
    on.exit(close(matlab))
  }

  R.matlab::setVariable(matlab, A=A, b=b, lambda=lambda, x0=initial)
  R.matlab::evaluate(matlab, 'x = solver_SLOPE(A, b, lambda, x0);')
  result = R.matlab::getVariable(matlab, 'x')
  as.vector(result$x)
}