#' Foreign Final Demand Domestic Value Added
#'
#' @name ffddva
#' @param x A Leontief decomposed Inter-Country Input Output table as created by decompr, which should be post multiplied with final demand (using the parameter: post="final_demand")
#' @param aggregate should dfddva be aggregated along source industries to a national sum?
#' @export
#' @import decompr
#' @examples
#' # load the decompr package
#' library(decompr)
#'
#' # load example data
#' data(leather)
#'
#' # create a leontief decomposed data set
#' l <- decomp(inter,
#'             final,
#'             countries,
#'             industries,
#'             out,
#'             method = "leontief",
#'             post = "final_demand")
#'
#'  # apply ffddva
#'  ffddva( l )

ffddva <- function ( x, aggregate=FALSE ) {

  # read attributes
  k      <- attr(x, "k")
  i      <- attr(x, "i")
  # rownam <- attr(x, "rownam")
  G  <- length(k)
  N  <- length(i)
  GN <- G*N

  # transform back to 2dim x 2dim matrix
  x <- matrix(x[,4], nrow=GN, byrow=TRUE)

  # remove exports to self
  x <- x - diagonals::fatdiag(diagonals::fatdiag( x, steps=G ), steps=G, nrow=GN, ncol=G )

  # sum across rows
  x <- rowSums(x)

  # create output format
  x <- data.frame(Source_Country = rep(k, each=N), Source_Industry = rep(i, times=G), ffddva = x)

  # aggregate or not
  if (aggregate) {
    f <- as.factor(x[,1])

    x <- tapply(x[,3], f, sum)

    x <- data.frame(Source_Country = row.names(x), ffddva = x)

    row.names(x) <- NULL
  }

  return(x)

}
