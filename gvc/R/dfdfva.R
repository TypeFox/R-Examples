#' Domestic Final Demand Foreign Value Added
#'
#' @name dfdfva
#' @param x A Leontief decomposed Inter-Country Input Output table as created by decompr, which should be post multiplied with final demand (using the parameter: post="final_demand")
#' @param aggregate should dfddva be aggregated along source industries to a national sum?
#' @export
#' @import decompr
#' @examples
#' # load the decompr package
#' library(decompr)
#'
#' # load the example data
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
#'  # apply dfdfva
#'  dfdfva( l )

dfdfva <- function ( x, aggregate=FALSE ) {

  # read attributes
  k      <- attr(x, "k")
  i      <- attr(x, "i")
  # rownam <- attr(x, "rownam")
  G <- length(k)
  N <- length(i)
  GN <- G*N

  # transform back to 2dim x 1dim matrix
  x <- matrix(x[,4], nrow=GN, ncol=G, byrow=TRUE)

  # remove everything except exports to self
  x <- x - diagonals::fatdiag(diagonals::fatdiag( x, steps=G ), steps=G, nrow=GN, ncol=G )

  # aggregate or not
  if (aggregate) {
    x <- colSums(x)

    x <- data.frame(Importing_Country = k , dfdfva = x)

    return(x)

  } else {

    # create factors for industries
    f <- gl(N, 1, GN)

    # create temporary matrix
    t <- matrix(0, nrow=N, ncol=G)

    for (j in 1:G) {
      t[,j] <- tapply(x[,j], f, sum)
    }

    x <- matrix(t, byrow=FALSE)

    x <- data.frame(Importing_Country = rep(k, each=N), Source_Industry = rep(i, times=G), dfdfva = x)

    return(x)

  }

}
