#' Extract spatial trend data
#' 
#' Extract spatial trend data from an object of class \code{likfit}.
#' 
#' @param x Object of class \code{likfit}.
#' 
#' @details 
#' 
#' \code{trend.terms} is similar to \code{\link[stats]{terms}}.
#' 
#' \code{trend.matrix} is similar to \code{\link[stats]{model.frame}}.
#' 
#' @seealso \code{\link[geoR]{likfit}}
#' 
#' @author Alessandro Samuel-Rosa \email{alessandrosamuelrosa@@gmail.com}
#' @aliases trend.terms trend.matrix
# FUNCTION - EXTRACT TREND TERMS ###############################################
# This function is similar to stats::terms
#' @export
#' @rdname trend
trend.terms <- 
  function (x) {
    cl <- class(x)
    if (all(cl == c("likGRF", "variomodel"))) {
      res <- all.vars(x$trend)
    }
    return (res)
  }
# FUNCTION - EXTRACT TREND MATRIX ##############################################
# This function is similar to stats::model.frame
#' @export
#' @rdname trend
trend.matrix <-
  function (x) {
    cl <- class(x)
    if (all(cl == c("likGRF", "variomodel"))) {
      res <- x$trend.matrix[, -1]
      colnames(res) <- trend.terms(x)
    }
    return (res)
  }
