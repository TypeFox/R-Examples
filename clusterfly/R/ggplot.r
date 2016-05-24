#' Static plot: Parallel coordinates.
#' Draw a parallel coordinates plot, facetted across clustering.
#'
#' This really only a proof of concept, a truly useful PCP
#' needs interaction, especially to move the variables around.
#'
#' @param cfly clusterfly object
#' @param index clustering to use
#' @param ... other arguments passed to \code{\link[ggplot2]{geom_line}}
#' @export
#' @keywords hplot
#' @examples
#' if (require("ggplot2")) {
#' o <- olive_example()
#' cfly_pcp(o, "kmeans")
#' cfly_pcp(o, "kmeans", alpha = 1/10)
#' cfly_pcp(o, "kmeans", alpha = 1/10) + coord_flip()
#' }
cfly_pcp <- function(cfly, index, ...) {
  stopifnot(require("ggplot2"))

  df <- data.frame(
    rescaler(cfly$df),
    .cluster = cfly$clusters[[index]],
    .id = 1:nrow(cfly$df))
  dfm <- melt(df, id = c(".cluster", ".id"))

  ggplot(dfm, aes_string(x = "variable", y = "value", group = ".id")) +
    geom_line(...) +
    facet_wrap(~ .cluster)
}

#' Static plot: Variable distribution.
#' Draw a density plot for each continuous variable, facetted across clustering.
#'
#' This allows you to quickly visualise how the cluster
#' vary in a univariate manner.  Currently, it is a bit
#' of a hack, because \code{\link[ggplot2]{ggplot}} does
#' not support plots with different scales, so the variables
#' are manually rescaled prior to plotting.
#'
#' This plot is inspired by Gaguin \url{http://www.rosuda.org/gaguin}.
#'
#' @param cfly clusterfly object
#' @param index clustering to use
#' @param scale scaling to use
#' @keywords hplot
#' @export
#' @examples
#' if (require("ggplot2")) {
#' o <- olive_example()
#' cfly_dist(o, "kmeans")
#' cfly_dist(o, "kmeans") + scale_y_continuous(limit=c(0, 2))
#' }
#' @importFrom reshape2 melt
cfly_dist <- function(cfly, index, scale="range") {
  stopifnot(require("ggplot2"))

  df <- cbind(cfly$df, .cluster=factor(cfly$clusters[[index]]))
  dfm <- melt(rescaler(df, scale), id=".cluster")

  ggplot(dfm, aes_string(x = "value")) +
    geom_density() +
    facet_grid(.cluster ~ variable)
}

#' Static plot: Fluctuation diagram.
#' Draw a fluctuation diagram comparing two clusterings.
#'
#' @param cfly clusterfly object
#' @param a first clustering, will be reordered to match \code{b} if clarify=TRUE
#' @param b second clustering
#' @param clarify use \code{\link{clarify}} to rearranged cluster indices?
#' @keywords hplot
#' @export
#' @importFrom plyr count
#' @examples
#' if (require("ggplot2")) {
#' o <- olive_example()
#' cfly_fluct(o, "kmeans", "Region")
#' cfly_fluct(o, "kmeans", "Region", clarify = FALSE)
#' }
cfly_fluct <- function(cfly, a, b, clarify=TRUE) {
  stopifnot(require("ggplot2"))

  ca <- cfly$clusters[[a]]
  cb <- cfly$clusters[[b]]
  if (clarify) ca <- clarify(ca, cb)

  counts <- count(data.frame(x = factor(ca), y = factor(cb)))
  nx <- length(levels(counts$x))
  ny <- length(levels(counts$y))

  counts$freq <- sqrt(counts$freq / max(counts$freq))

  ggplot(counts, aes_string(x = "x", y = "y", height = "freq", width = "freq")) +
    geom_tile(colour = "white") +
    scale_y_discrete(a) +
    scale_x_discrete(b) +
    theme(aspect.ratio = ny/nx)
}
