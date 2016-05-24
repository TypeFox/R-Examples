#' K-means on the Overall Mean, Row Margins or Column Margins
#' 
#' This function conducts k-means on the overall mean, the row margins or column margins of a set
#' of N matrices. These matrices are two-way slices of a three-dimensional array. 
#' 
#' @param data A three-way array representing the data.
#' @param margin An integer giving the single subscript of \code{data} over which the clustering 
#' will be applied. 
#' @param delta A four-element binary vector (logical or numeric) indicating which sum-to-zero 
#' constraints must be enforced.
#' @param nclust An integer giving the desired number of clusters. In case \code{type} specifies 
#' more than one method, \code{nclust} can be a vector containing the number of 
#' clusters to be determined for each type of cluster, and in the correct order as determined by
#' \code{type} (after matching the arguments). If \code{type} is of length greater than one and 
#' \code{nclust} is of length one, the behaviour is governed by \code{sep.nclust}.  
#' @param sep.nclust Logical indicating how nclust should be used across different \code{type}'s.
#' If \code{sep.nclust} is \code{TRUE}, \code{nclust} is recycled so that each \code{type} can
#' have a different number of clusters. If \code{sep.nclust} is \code{FALSE}, the same vector
#' \code{nclust} is used for all \code{type}'s.
#' @param type One of \code{"overall"}, \code{"rows"} or \code{"columns"} (or a unique abbreviation of 
#' one of these) indicating whether clustering should be done on row margins, column margins or
#' the overall means of the two-way slices respectively. If more than one opion are supplied, the
#' algorithm is run for all (unique) options supplied.
#' @param verbose Integer controlling the amount of information printed: 0 = no information, 
#' 1 = Information on random starts and progress, and 2 = information is printed after
#' each iteration for the interaction clustering.
#' @param \dots Additional arguments passed to \code{\link{kmeans}}.
#' @return A list containing a subset of the classes \code{row.kmeans}, \code{col.kmeans} and 
#' \code{ovl.kmeans} which are specific versions of class \code{kmeans}. In case \code{type} is a vector, a list
#' is returned containing the results for each of the (unique) elements of \code{type}, with the
#' same classes as before. See \code{\link{kmeans}} for an overview of the structure of these objects.
#' @seealso \code{\link{kmeans}}
#' @aliases ovl.kmeans row.kmeans col.kmeans
#' @export
orc.lsbclust <- function(data, margin = 3L, delta, nclust, sep.nclust = TRUE,  
                         type = NULL, verbose = 1, ...){
  
  ## Sanity checks and coercion
  delta <- as.numeric(delta)
  if (length(delta) != 4 || !all(delta %in% 0:1))  stop("Argument 'delta' supplied in an unrecognized format.")
  if (!is(data, "array") || length(dim(data)) != 3) stop("Data must be a three-way array.") 
  if (!all(margin %in% 1:3) || length(margin) != 1) stop("Argument 'margin' must be 1, 2 or 3.")
  if (all(delta[1:2] == 0)) return()
  
  ## Determine type from delta if not explicitly supplied
  if (is.null(type)) {
    ind.type <- rep(FALSE, 3L)
    ind.type[1] <- (delta[1] * delta[3] + delta[2] * delta[4] - delta[1] * delta[2]) == 1
    ind.type[2] <- delta[2] == 1
    ind.type[3] <- delta[1] == 1
    type <- c("overall", "rows", "columns")[ind.type]
    if(length(nclust) == 3) nclust <- nclust[ind.type]
  } else {
    type <- match.arg(tolower(type), choices = c("overall", "rows", "columns"), several.ok = TRUE)
    type <- unique(type)
  }
  
  ## Recurse if nclust is a vector
  if(length(nclust) > 1 && !sep.nclust) {
    out <- lapply(nclust, orc.lsbclust, data = data, margin = margin, delta = delta,
                  type = type, verbose = verbose, ...)
    class(out) <- "orc.kmeans"
    return(out)
  }
  
  ## Recurse if more than one 'type' is supplied
  if(length(type) > 1) {
    out <- mapply(orc.lsbclust, type = type, nclust = nclust, 
                  MoreArgs = list(data = data, margin = margin, delta = delta, verbose = verbose, ...), 
                  SIMPLIFY = FALSE)
    class(out) <- "orc.kmeans"
    return(out)
  }
  
  ## Return NULL when component is not to be fit according to delta
  if (delta[2] == 0 && type == "rows") return()
  if (delta[1] == 0 && type == "columns") return()
  if (delta[1] * delta[3] + delta[2] * delta[4] - delta[1] * delta[2] == 0 && type == "overall") return()
  
  ## Determine summarizing function for two-way slices
  sfun <- switch(type, rows = rowMeans, columns = colMeans, overall = mean)
  
  ## Determine post-summarizing centring functions
  colCentre <- function(x) sweep(x, MARGIN = 2L, STATS = colMeans(x), FUN = "-")
  postfun <- switch(type, rows = switch(delta[4L] + 1L, identity, colCentre), 
                         columns = switch(delta[3L] + 1L, identity, colCentre),
                         overall = identity)
  
  ## Calculate two-way data to be passed to kmeans()
  tdata <- apply(data, MARGIN = margin, FUN = sfun)
  tdata <- postfun(tdata)
  if(is.null(dim(tdata))) tdata <- matrix(tdata, ncol = 1)
  else tdata <- t(tdata)
  
  # Message if verbose
  if (verbose) {
    mess <- switch(type, rows = "K-means on row margins...", 
                   columns = "K-means on column margins...",
                   overall = "K-means on overall means...")
    cat(mess)
  }
  
  ## Try running k-means and process results
  out <- try(kmeans(tdata, centers = nclust, ...))
  
  ## Complete message
  if (verbose) {
    cat("\tDONE\n")
  }
  
  ## Update output unless error was encountered
  if (!inherits(out, "try-error")) {
    ## Degrees-of-freedom
    out$df <- switch(type, rows = delta[2] * (nclust * (ncol(tdata) - delta[4]) + nrow(tdata) * (nclust - 1)),
                     columns = delta[1] * (nclust * (ncol(tdata) - delta[3]) + nrow(tdata) * (nclust - 1)),
                     overall = (delta[1] * delta[3] + delta[2] * delta[4] - delta[1] * delta[2]) * (nclust + nrow(tdata) * (nclust - 1)))
    
#     ## Calinski-Harabasz statistic
#     N <- dim(data)[margin]
#     out$chstat <- out$betweenss * (N - nclust) / (out$tot.withinss * (nclust - 1))
    
    ## Reorder classes according to size
    ord <- order(out$size, decreasing = TRUE)
    out$size <- out$size[ord]
    out$centers <- out$centers[ord, , drop = FALSE]
    rownames(out$centers) <- seq_len(nclust)
    out$withinss <- out$withinss[ord]
    out$cluster <- plyr::mapvalues(x = out$cluster, from = ord, to = seq_along(ord))
    
    ## Set class and return
    outclass <- switch(type, rows = c("row.kmeans", "kmeans"), columns = c("col.kmeans", "kmeans"), 
                       overall = c("ovl.kmeans", "kmeans"))
    class(out) <- outclass
  }
  
  return(out)
}