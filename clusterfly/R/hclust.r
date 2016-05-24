# Need a new type of linking to make this work
# Brushing a node should highlight all nodes and leaves below it
# (investigate nested set representation for efficient storage)
#
# Can this be done using rggobi and the old linking code?
# Should it be added to ggobi as a new type of linking?
#    * If so, as the general case - defined by edges?
#    * Or for the particular nested set representation?


#' Visualisig hierarchical clustering.
#' This method supplements a data set with information needed to draw a
#' dendrogram
#'
#' Intermediate cluster nodes are added as needed, and positioned at the
#' centroid of the combined clusters.
#'
#' @param data data set
#' @param metric distance metric to use, see \code{\link{dist}} for list of
#'   possibilities
#' @param method cluster distance measure to use, see \code{\link{hclust}} for
#'   details
#' @return object of type, hierfly
#' @seealso \code{\link{cut.hierfly}}, \code{\link{ggobi.hierfly}}
#' @keywords cluster
#' @export
#' @examples
#' h <- hierfly(iris)
#' ggobi(h)
#' h <- hierfly(iris, method="single")
hierfly <- function(data, metric="euclidean", method="average") {
  cat <- sapply(data, is.factor)
  h <- hclust(dist(data[,!cat], metric), method)

  data$ORDER <- order(h$order)
  data$HEIGHT <- 0
  data$LEVEL <- 0
  data$POINTS <- 1

  for (i in 1:nrow(h$merge)) {
    newr <- combinerows(data[as.character(-h$merge[i,]),], cat)
    newr$HEIGHT <- h$height[i]
    newr$LEVEL <- i
    rownames(newr) <- as.character(-i)

    data <- rbind(data, newr)
  }

  data$node <- (as.numeric(rownames(data)) < 0) + 0

  structure(list(data=data, hclust=h), class="hierfly")
}

combinerows <- function(df, cat) {
  same <- function(x) if (length(unique(x)) == 1) x[1] else NA
  points <- df$POINTS

  cont <- as.data.frame(lapply(df[, !cat, drop=FALSE] * points, sum)) / sum(points)
  cat <- as.data.frame(lapply(df[, cat, drop=FALSE], same))

  df <- if (nrow(cont) > 0 && nrow(cat) > 0) {
    cbind(cont, cat)
  } else if (nrow(cont) > 0) {
    cont
  } else {
    cat
  }
  df$POINTS <- sum(points)
  df
}

#' @export
print.hierfly <- function(x, ...) {
  print(str(x))
}

#' Visualise hierarchical clustering with GGobi.
#' Displays both data and dendrogram in original high-d space.
#'
#' This adds four new variables to the original data set:
#'
#' \itemize{
#'    \item ORDER, the order in which the clusters are joined
#'    \item HEIGHT, the height of the branch, ie. the dissimilarity between the branches
#'    \item LEVEL, the level of the branch
#'    \item POINTS, the number of points in the branch
#' }
#'
#' Make sure to select "attach edge set (edges)" in the in the edges menu on the
#' plot window, when you create a new plot.
#'
#' A tour over the original variables will show how the clusters agglomerate
#' in space.   Plotting order vs height, level or points will give various
#' types of dendograms.  A correlation tour with height/level/points on the y
#' axis and the original variables on the x axis will show a mobile blowing
#' in the wind.
#'
#' @param data hierfly object to visualise in GGobi
#' @param ... ignored
#' @seealso \code{\link{cut.hierfly}}
#' @keywords cluster dynamic
#' @export
#' @examples
#' h <- hierfly(iris)
#' ggobi(h)
#' h <- hierfly(iris, method="single")
ggobi.hierfly <- function(data, ...) {
  h <- data$hclust
  data <- data$data

  g <- ggobi(data)
  d <- g[1]
  glyph_type(d) <- ifelse(data$node != 0, 1, 6)

  e <- data.frame(level=1:length(h$height), height=h$height)[rep(1:length(h$height), 2), ]
  rownames(e) <- paste("e", 1:nrow(e), sep="")

  g$edges <- e
  edges(g$edges) <- cbind(as.character(-h$merge), -rep(1:nrow(h$merge), 2))

  d <- displays(g)[[1]]
  edges(d) <- g[2]

  invisible(g)
}

#' Cut hierfly object into k clusters/colours.
#'
#' @param x hierfly object to colour
#' @param k number of clusters
#' @param g GGobi instance displaying x, will create new if not specified
#' @param ... ignored
#' @keywords cluster
#' @export
#' @examples
#' h <- hierfly(iris)
#' hfly <- ggobi(h)
#' cut(h, 2, hfly)
#' h <- hierfly(iris, method="ward")
#' g <- ggobi(h)
#' cut(h, 2, g)
cut.hierfly <- function(x, k=2, g=ggobi(x), ...) {
  d <- g[1]
  glyph_colour(d) <- c(cutree(x$hclust, k=k) + 1, rep(1, length(x$hclust$height)))
}
