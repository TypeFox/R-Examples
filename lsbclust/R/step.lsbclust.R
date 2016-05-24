#' Model Search for lsbclust
#' 
#' Fit \code{\link{lsbclust}} models for different numbers of clusters and/or different values of 
#' \code{delta}. The resulting output can be inspected through its \code{plot} method to facilitate 
#' model selection. Each component of the model is fitted separately.
#' 
#' @param data A three-way array representing the data.
#' @param margin An integer giving the single subscript of \code{data} over which the clustering 
#' will be applied. 
#' @param delta A four-element binary vector (logical or numeric) indicating which sum-to-zero 
#' constraints must be enforced.
#' @param nclust Either a vector giving the number of clusters which will be applied to each element
#' of the model, that is to (a subset of) the overall mean, row margins, column margins and 
#' interactions. If it is a list, arguments are matched by the names \code{"overall"}, \code{"rows"}
#' \code{"columns"} and \code{"interactions"}. If the list does not have names, the components are 
#' extracted in the aforementioned order.
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
#' @param verbose The number of iterations after which information on progress is provided 
#' (passed to \code{\link{int.lsbclust}}).
#' @param type One of \code{"rows"}, \code{"columns"} or \code{"overall"} (or a unique abbreviation of 
#' one of these) indicating whether clustering should be done on row margins, column margins or
#' the overall means of the two-way slices respectively. If more than one opion are supplied, the
#' algorithm is run for all (unique) options supplied (passed to \code{\link{orc.lsbclust}}). This
#' is an optional argument.
#' @param \dots Additional arguments passed to \code{\link{kmeans}}.
#' @export
#' @examples
#' m <- step.lsbclust(data = dcars, margin = 3, delta = c(1, 0, 1, 0), nclust = 4:5, 
#'                      ndim = 2, fixed = "columns", nstart = 1, nstart.kmeans = 100, 
#'                      parallelize = FALSE)
#'                      
#' ## For a list of all deltas                     
#' delta <- expand.grid(replicate(4, c(0,1), simplify = FALSE))
#' delta <- with(delta, delta[!(Var1 == 0 & Var3 == 1), ])
#' delta <- with(delta, delta[!(Var2 == 0 & Var4 == 1),])
#' delta <- delta[-4,]
#' delta <- as.list(as.data.frame(t(delta)))
#' m2 <- step.lsbclust(data = dcars, margin = 3, delta = delta, nclust = 4:5, 
#'                      ndim = 2, fixed = "columns", nstart = 1, nstart.kmeans = 100, 
#'                      parallelize = FALSE)
step.lsbclust <- function(data, margin = 3L, delta = c(1, 1, 1, 1), nclust, ndim = 2,
                          fixed = c("none", "rows", "columns"), nstart = 20, starts = NULL, 
                          nstart.kmeans = 500, alpha = 0.5, parallelize = FALSE, maxit = 100, 
                          verbose = -1, type = NULL,  ...) {
  
  ## Start timing
  time0 <- proc.time()[3]
  
  ## Recurse of delta is a list
  if (is.list(delta)) {
    return(lapply(delta, step.lsbclust, data = data, margin = margin, nclust = nclust, ndim = ndim, 
                  fixed = fixed, nstart = nstart, starts = starts, nstart.kmeans = nstart.kmeans, 
                  alpha = alpha, parallelize = parallelize, maxit = maxit, verbose = verbose, 
                  type = type, ...))
  }
  
  ## Capture call
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
  
  ## If is nclust is vector, create a named list of 4 elements out of it
  if (!is.list(nclust)) {
    if (is.numeric(nclust)) {
      nclust <- replicate(4, nclust, simplify = FALSE)
      names(nclust) <- c("overall", "rows", "columns", "interactions")
    }
  } else {
    ## Otherwise, check length and names
#     if (length(nclust) != 4) stop("nclust is not of length 4.")
    if (!all(names(nclust) %in% c("overall", "rows", "columns", "interactions")) || is.null(names(nclust))) 
      names(nclust) <- c("overall", "rows", "columns", "interactions")
  }
  
  ## Setup up return values
  res <- vector(length = 4, mode = "list")
  names(res) <- names(nclust)
  
  ## Do the overall means
  if (!is.null(nclust$overall)) 
    res$overall <- orc.lsbclust(data = data, margin = margin, delta = delta, nclust = nclust$overall,
                              type = "overall", nstart = nstart.kmeans, sep.nclust = FALSE, ...)
  
  ## Do the rows
  if (!is.null(nclust$rows)) 
    res$rows <- orc.lsbclust(data = data, margin = margin, delta = delta, nclust = nclust$rows,
                              type = "rows", nstart = nstart.kmeans, sep.nclust = FALSE, ...)
  
  ## Do the overall means
  if (!is.null(nclust$columns)) 
    res$columns <- orc.lsbclust(data = data, margin = margin, delta = delta, nclust = nclust$columns,
                              type = "columns", nstart = nstart.kmeans, sep.nclust = FALSE, ...)
  
  ## Do the interactions
  if (!is.null(nclust$interactions))
      res$interactions <- int.lsbclust(data = data, margin = margin, delta = delta, nclust = nclust$interactions, 
                                   ndim = ndim, fixed = fixed, nstart = nstart, starts = starts, 
                                   alpha = alpha, parallelize = parallelize, maxit = maxit, 
                                   verbose = verbose, method = NULL)
  
  ## Extract losses
  loss <- vector(length = 4, mode = "list")  
  names(loss) <- names(nclust)
  if (delta[1] * delta[3] + delta[2] * delta[4] - delta[1] * delta[2]) {
    loss$overall <- J * K * sapply(res$overall, "[[", "tot.withinss")
  }
  if (delta[2]) loss$rows <- K * sapply(res$rows, "[[", "tot.withinss")
  if (delta[1]) loss$columns <- J * sapply(res$columns, "[[", "tot.withinss")
  loss$interactions <- sapply(res$interactions, "[[", "minloss") * sapply(res$interactions, "[[", "maxloss")
  
  ## Extract df's
  df <- vector(length = 4, mode = "list")  
  names(df) <- names(nclust)
  if (delta[1] * delta[3] + delta[2] * delta[4] - delta[1] * delta[2]) {
    df$overall <- sapply(res$overall, "[[", "df")
  }
  if (delta[2]) df$rows <- sapply(res$rows, "[[", "df")
  if (delta[1]) df$columns <- sapply(res$columns, "[[", "df")
  df$interactions <- sapply(res$interactions, "[[", "df") 

  ## Construct all combinations of losses
  ind <- (1:4)[!sapply(loss, is.null)]
  losscomb <- expand.grid(loss[ind])
  losscomb$total <- rowSums(losscomb)

  ## Construct corresponding df
  dfcomb <- expand.grid(df[ind])
  dfcomb$total <- rowSums(dfcomb)
  
  ## Corresponding indices for nclust
  nclustcomb <- expand.grid(nclust[ind])
  
  ## Return
  out <- list(loss = losscomb, df = dfcomb, nclust = nclustcomb, 
              ind.loss = loss[ind], ind.df = df[ind], ind.nclust = nclust[ind], time = proc.time()[3] - time0)
  class(out) <- "step.lsbclust"
  return(out)
}