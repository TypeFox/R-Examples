#' Leontief Decomposition
#'
#' @param x an object of class decompr
#' @param post post-multiply the Leontief inverse with something, the default is exports
#' @param long transform the output data into a long (tidy) data set or not, default it TRUE.
#' @return a data frame containing the square matrix and labelled column and rows
#' @author Bastiaan Quast
#' @references Wang, Zhi, Shang-Jin Wei, and Kunfu Zhu.
#' Quantifying international production sharing at the bilateral and sector levels.
#' No. w19677. National Bureau of Economic Research, 2013.
#' @export
#' @examples
#' # load example data
#' data(leather)
#' 
#' # create intermediate object (class decompr)
#' decompr_object <- load_tables_vectors(inter,
#'                                       final,
#'                                       countries,
#'                                       industries,
#'                                       out        )
#'
#' # run the Leontief decomposition on the decompr object
#' leontief(decompr_object )


leontief <- function( x, post = c("exports", "output", "final_demand", "none"), long=TRUE ) {

  post <- match.arg(post)

  # compute Leontief inverse
  out <- x$Vhat %*% x$B

  # post multiply
  if (post == "exports") {
    out <- out %*% x$Exp
  } else if (post == "output") {
    out <- out %*% diag(x$X)
  } else if (post == "final_demand") {

    out <- out %*% x$Y

    # create output format for post="final_demand"
    out <- as.vector(t(out))
    out <- data.frame( rep(x$k, each=x$GN),
                       rep(x$i, times=x$G, each=x$G),
                       rep(x$k, times=x$GN),
                       out
                       )
    names(out) <- c("Source_Country", "Source_Industry", "Importing_Country", "Final_Demand")

    # create attributes
    attr(out, "k")      <- x$k
    attr(out, "i")      <- x$i
    attr(out, "decomposition") <- "leontief"

    return(out)
  }



  # structure output format
  if (long == TRUE) {

    out <- as.vector(t(out))
    out <- data.frame( rep(x$k,                  each = x$GN*x$N ),
                       rep(x$i, times = x$G,     each = x$GN),
                       rep(x$k, times = x$GN,    each = x$N),
                       rep(x$i, times = x$GN*x$G ),
                       out)
    names(out) <- c("Source_Country", "Source_Industry", "Using_Country", "Using_Industry", "FVAX")

    # set long attribute to TRUE
    attr(out, "long") <- TRUE

  } else {

    # add row and column names
    out <- as.data.frame(out)
    if (post != "final_demand") names(out) <- x$rownam
    row.names(out) <- x$rownam

    # set long attribute to FALSE
    attr(out, "long") <- FALSE

  }



  # create attributes
  attr(out, "k")      <- x$k
  attr(out, "i")      <- x$i
  attr(out, "decomposition") <- "leontief"
  # attr(out, "rownam") <- x$rownam

  # return result
  return( out )

}
