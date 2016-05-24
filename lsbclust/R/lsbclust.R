#' Least-squares Bilinear Clustering of Three-way Data
#' 
#' This function clusters along one way of a three-way array (as specified by \code{margin}) while
#' decomposing along the other two dimensions. Four types of clusterings are allowed based on the
#' respective two-way slices of the array: on the overall means, row margins, column margins and the 
#' interactions between rows and columns. Which clusterings can be fit is determined by the vector
#' \code{delta}, with four binary elements. All orthogonal models are fitted. 
#' The nonorthogonal case \code{delta = (1, 1, 0, 0)} returns an error. See the reference for further details.
#' 
#' @param data A three-way array representing the data.
#' @param margin An integer giving the single subscript of \code{data} over which the clustering 
#' will be applied. 
#' @param delta A four-element binary vector (logical or numeric) indicating which sum-to-zero 
#' constraints must be enforced.
#' @param nclust A vector of length four giving the number of clusters for the overall mean, the row
#' margins, the column margins and the interactions (in that order) respectively. Alternatively, a
#' vector of length one, in which case all components will have the same number of clusters.
#' @param ndim The required rank for the approximation of the interactions (a scalar).
#' @param fixed One of \code{"none"}, \code{"rows"} or \code{"columns"} indicating whether to fix neither
#' sets of coordinates, or whether to fix the row or column coordinates across clusters respectively.
#' If a vector is supplied, only the first element will be used (passed to \code{\link{int.lsbclust}}).
#' @param nstart The number of random starts to use for the interaction clustering.
#' @param nstart.kmeans The number of random starts to use in \code{\link{kmeans}}.
#' @param starts A list containing starting configurations for the cluster membership vector. If not
#' supplied, random initializations will be generated (passed to \code{\link{int.lsbclust}}).
#' @param alpha Numeric value in [0, 1] which determines how the singular values are distributed
#' between rows and columns (passed to \code{\link{int.lsbclust}}).
#' @param parallelize Logical indicating whether to parallelize over different starts or not 
#' (passed to \code{\link{int.lsbclust}}).
#' @param maxit The maximum number of iterations allowed in the interaction clustering.
#' @param verbose Integer controlling the amount of information printed: 0 = no information, 
#' 1 = Information on random starts and progress, and 2 = information is printed after
#' each iteration for the interaction clustering.
#' @param method The method for calculating cluster agreement across random starts, passed on
#' to \code{\link{cl_agreement}} (passed to \code{\link{int.lsbclust}}).
#' @param sep.nclust Logical indicating how nclust should be used across different \code{type}'s.
#' If \code{sep.nclust} is \code{TRUE}, \code{nclust} is recycled so that each \code{type} can
#' have a different number of clusters. If \code{sep.nclust} is \code{FALSE}, the same vector
#' \code{nclust} is used for all \code{type}'s.
#' @param type One of \code{"rows"}, \code{"columns"} or \code{"overall"} (or a unique abbreviation of 
#' one of these) indicating whether clustering should be done on row margins, column margins or
#' the overall means of the two-way slices respectively. If more than one opion are supplied, the
#' algorithm is run for all (unique) options supplied (passed to \code{\link{orc.lsbclust}}). This
#' is an optional argument.
#' @param \dots Additional arguments passed to \code{\link{kmeans}}.
#' @return Returns an object of S3 class \code{lsbclust} which has slots:
#'    \item{\code{overall}}{Object of class \code{ovl.kmeans} for the overall means clustering}
#'    \item{\code{rows}}{Object of class \code{row.kmeans} for the row means clustering}
#'    \item{\code{columns}}{Object of class \code{col.kmeans} for the column means clustering}
#'    \item{\code{interactions}}{Object of class \code{int.lsbclust} for the interaction clustering}
#'    \item{\code{call}}{The function call used to create the object}
#'    \item{\code{delta}}{The value of \code{delta} in the fit}
#'    \item{\code{df}}{Breakdown of the degrees-of-freedom across the different subproblems}
#'    \item{\code{loss}}{Breakdown of the loss across subproblems}
#'    \item{\code{time}}{Time taken in seconds to calculate the solution}
#'    \item{\code{cluster}}{Matrix of cluster membership per observation for all cluster types}
#' @seealso \code{\link{int.lsbclust}}, \code{\link{orc.lsbclust}}
#' @export
#' @references Schoonees, P.C., Groenen, P.J.F., Van de Velden, M. Least-squares Bilinear Clustering
#' of Three-way Data. Econometric Institute Report, EI2014-23.
#' @importFrom graphics plot
#' @importFrom methods is
#' @import stats
lsbclust <- function(data, margin = 3L, delta = c(1, 1, 1, 1), nclust, ndim = 2,
                     fixed = c("none", "rows", "columns"), nstart = 20, starts = NULL, 
                     nstart.kmeans = 500, alpha = 0.5, parallelize = FALSE, maxit = 100, verbose = 1, 
                     method = "diag", type = NULL, sep.nclust = TRUE, ...) {
  
  ## Capture call, start time
  time0 <- proc.time()[3]
  cll <- match.call()
  
  ## Preliminaries
  J <- dim(data)[-margin][1]
  K <- dim(data)[-margin][2]
  
  ## Select correct option for fixed (only one)
  fixed <- match.arg(fixed)
  
  ## Sanity checks
  delta <- as.numeric(delta)
  if (length(delta) != 4 || !all(delta %in% 0:1))  stop("Argument 'delta' supplied in an unrecognized format.")
  if (all(delta == c(1, 1, 0, 0))) stop("This model is non-orthogonal")
  if (!is(data, "array") || length(dim(data)) != 3) stop("Data must be a three-way array.") 
  if (!all(margin %in% 1:3) || length(margin) != 1) stop("Argument 'margin' must be 1, 2 or 3.")
  
  ## Check also ... arguments against stats::kmeans()
  args.dots <- list(...)
  nms.dots <- names(args.dots)
  kmeans.args <- names(formals(stats::kmeans))
  check.dots <- nms.dots %in% kmeans.args
  if (any(!check.dots)) {
    stop("Unknown argument(s): ", nms.dots[!check.dots])    
  }
  
  ## Make sure there are dimnames for ways not clustered over
  dnms <- dimnames(data)
  if (is.null(dnms)) namenull <- rep(TRUE, 3L)
  else namenull <- sapply(dimnames(data), is.null)
  if (any(namenull[-margin])) {
    dims <- !(seq_len(3) == margin)
    dimnames(data)[dims & namenull] <- list(paste0("Row", seq_len(J)), paste0("Col", seq_len(K)))[namenull[-margin]]
  }
  
  ## If nclust is of length one, expand it to length 4
  if (length(nclust) == 1) nclust <- rep(nclust, 4)
  if (length(nclust) != 4) stop("nclust should be either of length 1 or 4.")
  
  ## Other components
  orc <- orc.lsbclust(data = data, margin = margin, delta = delta, nclust = nclust[-4L], 
                      sep.nclust = sep.nclust, type = type, nstart = nstart.kmeans, 
                      verbose = verbose, ...)
  if (all(delta == c(1, 0, 0, 0)) || all(delta == c(0, 1, 0, 0)) || 
        all(delta == c(0, 1, 1, 0)) || all(delta == c(1, 0, 0, 1))) {
    orc <- list(orc)
    names(orc) <- ifelse(delta[1], "columns", "rows")
  }
  ind.orc <- c("overall", "rows", "columns") %in% names(orc)
  
  ## Interactions
  if (verbose) {
    cat(paste0("Interaction clustering (", nstart, " starts)..."))
  }
  int <- int.lsbclust(data = data, margin = margin, delta = delta, nclust = nclust[4L], ndim = ndim,
                      fixed = fixed, nstart = nstart, starts = starts, alpha = alpha, 
                      parallelize = parallelize, maxit = maxit, verbose = verbose, method = method)
  if (verbose) {
    cat("\tDONE\n")
  }
  
  ## Total degrees-of-freedom
  df <- unlist(c(sapply(orc, '[[', 'df'), interactions = int$df))
  
  ## Total loss
  loss <- unlist(c(mapply(function(x, y) x$tot.withinss * y, orc, c(J*K, K, J)[ind.orc]), 
            interactions = int$minloss * int$maxloss))
  
  ## Process output
  out <- list(overall = orc$overall, rows = orc$rows, columns = orc$columns, interactions = int, 
              call = cll, delta = delta, df = df, loss = loss, time = proc.time()[3] - time0)
  
  ## Collect all clusterings
  clust <- do.call(cbind, lapply(out[1:4], "[[", "cluster"))
  out$cluster <- clust
  
  class(out) <- "lsbclust"
  return(out)
}