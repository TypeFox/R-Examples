#' @title Improved embed() function
#' @description 
#' Improved version of the \code{\link[stats]{embed}} function in 
#' the \code{stats} package.
#' First it allows embeddings in past and future observation space (backward and
#' forward shift). Secondly, it adds 'NA' to the beginning (or end) of the embedding
#' matrix, depending on the dimension of the embedding. Optionally, they can be removed.
#'
#' @param x a numeric vector, matrix, or time series.
#' @param max.lag a scalar representing the embedding dimension in past or future.
#' Note that contrary to 'dimension = 1' in \code{\link[stats]{embed}}, here 
#' 'max.lag = 1' will return a 2 column matrix (0 lag, 1 lag), 
#' and not just a 1 column matrix.  Similarly, for negative shift; e.g., 'max.lag = -2' returns
#' 3 column matrix with (0 lag, -1 lag, -2 lag).
#' @param na.omit logical; if TRUE, it removes NA values automatically from embedded matrix.
#' @export
#' @seealso \code{\link[stats]{embed}}
#' @examples
#' data(nottem)
#' aa <- embed2(nottem, 12)

embed2 <- function(x, max.lag = 1, na.omit = FALSE) {
  
  if (max.lag == 0) {
    embedded.space <- cbind(x)
  } else {
    if (max.lag < 0) {
      future <- TRUE
      if (is.null(dim(x))) {
        x <- rev(x)
      } else {
        x <- apply(x, 2, rev)
      }
      max.lag <- -max.lag
    } else {
      future <- FALSE
    }
    
    if (na.omit) {
      embedded.space <- embed(x, max.lag + 1)
    } else {
      if (is.null(dim(x))) {
        embedded.space <- embed(c(rep(NA, max.lag), x), max.lag + 1)
      } else {
        embedded.space <- embed(rbind(matrix(rep(NA, max.lag * ncol(x)),
                                             ncol = ncol(x)),
                                      x),
                                max.lag + 1)
      }
    }  
      
    if (future) {
      embedded.space <- apply(embedded.space, 2, rev)
      #colnames(embedded.space) <- paste0("fwd.", 0:max.lag)
    } else {
      #colnames(embedded.space) <- paste0("back.", 0:max.lag)
    }
  }
  return(embedded.space)
} 
