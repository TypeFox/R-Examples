# We compute the solution path of the fused lasso problem:
#
# \hat{\beta}(\lambda) =
# \argmin_\beta \|y - X \beta|_2^2 + \lambda\|D \beta\|_1,
#
# where D is the incidence matrix of a given graph, and X is a full
# column rank predictor matrix. The solution is piecewise constant
# over the graph.
#
# Note: if a matrix D is passed, but no graph, there is no error
# checking whatsoever; i.e. we don't check that D is actually an
# incidence matrix. We leave it up to the user to respect this.
#
# Note: we also do not check that X has full column rank.

fusedlasso <- function(y, X, D, graph, gamma=0, approx=FALSE,
                       maxsteps=2000, minlam=0, rtol=1e-7, btol=1e-7,
                       eps=1e-4, verbose=FALSE) {

  if (missing(y)) stop("y is missing.")
  if (!is.numeric(y)) stop("y must be numeric.")
  if (length(y) == 0) stop("There must be at least one data point [must have length(y) > 1].")

  if (missing(X)) X = NULL
  if (!is.null(X) && !is.matrix(X)) stop("X must be a matrix.")
  # Right now there is no check for X having full
  # column rank; this is for efficiency, we hence
  # never have to compute its pseudoinverse

  if (missing(D)) D = NULL
  if (missing(graph)) graph = NULL

  if (is.null(D) && is.null(graph)) {
    stop("Must pass either incidence matrix D or graph object.")
  }
  if (!is.null(D) && !is.matrix(D) && c(attributes(class(D))$package,"")[[1]] != "Matrix") {
    stop("D must be a matrix or a Matrix (from the Matrix package).")
  }
  if (!is.null(graph) && !any(class(graph)=="igraph")) {
    stop("Invalid graph object.")
  }
  if (!is.null(graph) && is.directed(graph)) {
    warning("Only undirected graphs are supported. Converting to undirected form.")
    graph = as.undirected(graph)
  }
  if (!is.null(D) && !is.null(graph)) {
    if (any(graph.laplacian(graph) != t(D)%*%D)) {
      stop("The incidence matrix D and the graph object are incompatable; specify only one.")
    }
  }

  # For simplicity
  y = as.numeric(y)

  # Build D if we have a graph
  if (is.null(D) && !is.null(graph)) {
    D = getDgSparse(graph)
  }

  # Make D a sparse matrix if it's not already
  if (c(attributes(class(D))$package,"")[[1]] != "Matrix") {
    D = Matrix(D,sparse=TRUE)
  }

  # Check dimensions
  if (is.null(X) && length(y)!=ncol(D)) stop("Dimensions don't match [length(y) != ncol(D)].")
  if (!is.null(X) && length(y)!=nrow(X)) stop("Dimensions don't match [length(y) != nrow(X)].")
  if (!is.null(X) && ncol(D)!=ncol(X)) stop("Dimensions don't match [ncol(X) != ncol(D)].")

  # Check gamma
  if (!is.numeric(gamma) || gamma<0) {
    stop("gamma must be between nonnegative.")
  }

  # Figure out whether we are applying a ridge penalty
  if (!is.null(X) && nrow(X) >= ncol(X)) eps=0
  else if (!is.null(X)) warning(sprintf(paste("Adding a small ridge penalty (multiplier %g),",
                                              "because X has more columns than rows."),eps))

  if (gamma==0) {
    if (is.null(X)) out = dualpathFused(y,D,approx,maxsteps,minlam,rtol,btol,verbose)
    else out = dualpathFusedX(y,X,D,approx,maxsteps,minlam,rtol,btol,eps,verbose)
  }
  else {
    # Construct the proper D matrix
    D0 = D
    D = rBind(D0,bandSparse(ncol(D0),k=0,diagonals=list(rep(gamma,ncol(D0)))))

    if (is.null(X)) out = dualpathFusedL1(y,D,D0,gamma,approx,maxsteps,minlam,rtol,btol,verbose)
    else out = dualpathFusedL1X(y,X,D,D0,gamma,approx,maxsteps,minlam,rtol,btol,eps,verbose)
  }

  out$call = match.call()

  class(out) = c("fusedlasso", "genlasso", "list")
  return(out)
}



