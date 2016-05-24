#' Rewrite a list of clauses to a string containing a Boolean formula in disjunctive normal form
#' 
#' @param dnf list of clauses
#' @return string containing a Boolean formula in disjunctive normal form
format_dnf <- function(dnf) {
  return(Reduce(function(x, y) {paste(x, y, sep = " + ")}, dnf))
}

#' Plot fuzzy set score of two sets against each other
#'
#' @param x Formula that describes the fuzzy set to plot along the x axis
#' @param y Formula that describes the fuzzy set to plot along the y axis
#' @param data Data set of basic fuzzy set scores
#' @param labels flag whether to label individual points with the case names
#' @param main.diagonal flag whether to plot the main diagonal
#' @param anti.diagonal flag whether to plot the anti diagonal
#' @return the \code{ggplot} plot object
#' 
#' @examples
#' require(QCAGUI)
#' data(d.urban)
#' xyplot("MLC", "WSR", d.urban)
#' 
#' @export
xyplot <- function(x, y, data, labels=FALSE, main.diagonal=TRUE, anti.diagonal=FALSE) {
 
  if (is.character(x)) {
    xname <- format_dnf(x)
  } else {
    xname <- ""
  }
  if (is.character(y)) {
    yname <- format_dnf(y)
  } else {
    yname <- ""
  }
  
  x <- evaluate_dnf(data, x)
  y <- evaluate_dnf(data, y)
  if (labels) {
    casenames <- rownames(data)
  }
  pl <- ggplot2::qplot(x,y) + ggplot2::scale_x_continuous(limits = c(0, 1), name = xname) +
                              ggplot2::scale_y_continuous(limits = c(0, 1), name = yname)
  if (labels) {
    pl <- directlabels::direct.label(pl + ggplot2::geom_point(ggplot2::aes(colour=casenames)))   
  }
  if (main.diagonal) {
    pl <- pl + ggplot2::geom_segment(ggplot2::aes(x = 0, y = 0, xend = 1, yend = 1))
  }
  if (anti.diagonal) {
    pl <- pl + ggplot2::geom_segment(ggplot2::aes(x = 0, y = 1, xend = 1, yend = 0))
  }
  return(pl)
}
  
#' Plot the fuzzy set scores of the solution and the outcome against each other
#'
#' @param x an object of class \code{qca} as returned by \code{\link{eqmcc}} of the package \code{QCA}
#' @param ... further arguments passed on to \code{\link{xyplot}}
#' @return the \code{ggplot} plot object
#' 
#' @examples
#' require(QCAGUI)
#' data(d.urban)
#' solution <- eqmcc(d.urban, outcome="RT", conditions=c("MLC", "FRB", "CP", "WSR"))
#' plot(solution)
#' 
#' @export
#' @S3method plot qca
plot.qca <- function(x, ...) {
  xyplot(toupper(x$tt$options$outcome), x$solution[[1]], x$tt$initial.data, ...)
}