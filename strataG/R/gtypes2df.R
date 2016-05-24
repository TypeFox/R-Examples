#' @title Convert \code{gtypes} To \code{data.frame}
#' @description Convert a \linkS4class{gtypes} object to a \code{data.frame}.
#'   
#' @param g a \linkS4class{gtypes} object.
#' @param one.col logical. If \code{TRUE}, then result has one column per 
#'   locus.
#' @param sep character to use to separate alleles if \code{one.col} is 
#'   \code{TRUE}.
#' @param ... additional arguments ot be passed to or from methods.
#'   
#' @return A \code{data.frame} with one row per sample.
#' 
#' @author Eric Archer \email{eric.archer@@noa.gov}
#' 
#' @seealso \link{df2gtypes}, \link{as.matrix.gtypes}
#' 
#' @export
#' 
gtypes2df <- function(g, one.col = TRUE, sep = "/", ...) {
  gen.df <- data.frame(as.matrix(g, one.col = one.col, sep = sep, ...), stringsAsFactors = FALSE)
  cbind(id = rownames(gen.df), strata = strata(g)[rownames(gen.df)], gen.df,
        stringsAsFactors = FALSE)
}