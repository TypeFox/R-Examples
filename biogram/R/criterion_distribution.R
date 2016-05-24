#' criterion_distribution class
#'
#' A result of \code{\link{distr_crit}} function.
#'
#' @details An object of class \code{criterion_distribution} is a numeric matrix. 
#' @section Data:
#' \describe{
#'   \item{1st column:}{possible values of criterion.}
#'   \item{2nd column:}{probability density function.}
#'   \item{3rd column:}{cumulative distribution function.}
#' }
#' @section Attributes:
#' \describe{
#'   \item{plot_data}{A matrix with values of the criterion and their probabilities.}
#'   \item{nice_name}{'Nice' name of the criterion.}
#' }
#' @name criterion_distribution
#' @docType class
NULL

create_criterion_distribution <- function(criterion, pdf, range, unsort_criterion,
                                          unsort_prob, nice_name) {
  dist <- cbind(criterion, 
                pdf, 
                1 - rev(cumsum(rev(pdf))))
  colnames(dist) <- c("criterion", "pdf", "cdf")

  attr(dist, "plot_data") <- matrix(c(unsort_criterion, unsort_prob), ncol = 2,
                                    dimnames = list(range, c("unsort_criterion",
                                                             "unsort_prob")))
  attr(dist, "nice_name") <- nice_name
  class(dist) <- "criterion_distribution"
  dist
}

#' Plot criterion distribution
#'
#' Plots results of \code{\link{distr_crit}} function.
#'
#' @param x object of class \code{\link{criterion_distribution}}.
#' @param ... further arguments passed to \code{\link[graphics]{plot}}.
#' @return nothing.
#' @export
#' @examples
#' target_feature <- create_feature_target(10, 375, 15, 600) 
#' example_result <- distr_crit(target = target_feature[,1], 
#'                              feature = target_feature[,2])
#' plot(example_result)
#' 
#' #a ggplot2 plot
#' library(ggplot2)
#' ggplot_distr <- function(x) {
#' b <- data.frame(cbind(x=as.numeric(rownames(attr(x, "plot_data"))), 
#'                       attr(x, "plot_data")))
#' d1 <- cbind(b[,c(1,2)], attr(x, "nice_name"))
#' d2 <- cbind(b[,c(1,3)], "Probability")
#' colnames(d1) <- c("x", "y", "panel")
#' colnames(d2) <- c("x", "y", "panel")
#' d <- rbind(d1, d2)
#' p <- ggplot(data = d, mapping = aes(x = x, y = y)) + 
#'   facet_grid(panel~., scale="free") + 
#'   geom_freqpoly(data= d2, aes(color=y), stat = "identity") + 
#'   scale_fill_brewer(palette = "Set1") + 
#'   geom_point(data=d1, aes(size=y), stat = "identity") + 
#'   guides(color = "none") + 
#'   guides(size = "none") + 
#'   xlab("Number of cases with feature=1 and target=1") + ylab("")
#' p
#' }
#' ggplot_distr(example_result)
#' 
plot.criterion_distribution <- function(x, ...) {
  old_par <- par(c("mar", "fig", "oma"))
  par(mar = c(5,4,4,5) + 0.1)
  plot(as.numeric(rownames(attr(x, "plot_data"))), 
       attr(x, "plot_data")[, "unsort_criterion"], 
       col="red", 
       xlab = "Number of cases with feature=1 and target=1",
       ylab = attr(x, "nice_name"))
  par(new = TRUE)
  plot(as.numeric(rownames(attr(x, "plot_data"))), 
       attr(x, "plot_data")[, "unsort_prob"], type = "l", 
       col = "green", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
  axis(4)
  mtext("density",side = 4,line = 3)
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), 
      new = TRUE)
  plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
  legend("top", legend = c(attr(x, "nice_name"), "Probability"), xpd = TRUE, 
         horiz = TRUE, fill = c("red", "green"), bty = "n", cex = 1)
  par(mar = old_par[["mar"]], fig = old_par[["fig"]], oma = old_par[["oma"]])
}

