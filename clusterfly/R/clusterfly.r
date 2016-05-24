#' Creates a convenient data structure for dealing with a dataset and a number
#' of alternative clusterings.
#'
#' Once you have created a clusterfly object, you can add
#' clusterings to it with \code{\link{cfly_cluster}}, and
#' visualise then in GGobi with \code{\link{cfly_show}} and
#' \code{\link{cfly_animate}}. Static graphics are also
#' available: \code{\link{cfly_pcp}} will produce a parallel
#' coordinates plot, \code{\link{cfly_dist}} will show
#' the distribution of each variable in each cluster, and
#' \code{\link{cfly_fluct}} compares two clusterings with a
#' fluctuation diagram.
#'
#' If you want to standardise the cluster labelling to one
#' group, look at \code{\link{clarify}} and \code{\link{cfly_clarify}}
#'
#' @param df data frame to be clustered
#' @param extra extra variables to be included in output, but not clustered
#' @param rescale rescale, if true each variable will be scaled to have mean 0
#'   and variance 1.
#' @seealso vignette("introduction")
#' @export
#' @aliases clusterfly package-clusterfly
#' @import rggobi
#' @keywords dynamic
#' @examples
#' ol <- olive_example()
#'
#' if (interactive()) {
#' ggobi(ol)
#' cfly_show(ol, "k4-1")
#' cfly_animate(ol, max = 5)
#' close(ol)
#' }
clusterfly <- function(df, extra = NULL, rescale=TRUE) {
  if (rescale) df <- rescaler(df)

  g <- NULL
  open_ggobi <- function() {
    if (is.null(g)) {
      clusters <- do.call("cbind", compact(list(df, extra)))
      g <<- ggobi(clusters)
    }
    invisible(g)
  }
  close_ggobi <- function() {
    if (is.null(g)) return()
    close(g)
    g <<- NULL
  }

  structure(list(
    df = df,
    extra = extra,
    clusters = list(),
    ggobi = open_ggobi,
    close = close_ggobi
  ), class="clusterfly")
}


#' Show in ggobi.
#' Opens an instance ggobi for this dataset (if not already open), and colours
#' the points according the cluster assignment.
#'
#' @param cf clusterfly object
#' @param idx clustering to display
#' @param hulls add convex hull? see \code{\link{addhull}} for details
#' @keywords dynamic
#' @export
#' @examples
#' o <- olive_example()
#' cfly_show(o, 1)
#' cfly_show(o, "Region")
#' if (!interactive()) close(o)
cfly_show <- function(cf, idx = "true", hulls = FALSE) {
  g <- cf$ggobi()[1]
  cl <- cf$clusters[[idx]]
  glyph_colour(g) <- cl
  if (hulls) {
    addhull(g[1], g, cl)
    glyph_colour(g['hulls']) <- g['hulls']$id
  }
}

#' @export
ggobi.clusterfly <- function(data, ...) data$ggobi()
#' @export
close.clusterfly <- function(con, ...) con$close()

#' @export
"[[<-.clusterfly" <- function(x, i, value) {
  x$clusters[[i]] <- value
  x
}

#' @export
print.clusterfly <- function(x, ...) {
  cat("Data:     ", paste(names(x$df), collapse=", "), "  [", nrow(x$df), "x", ncol(x$df), "]\n", sep="")
  cat("Extra:    ", paste(names(x$extra), collapse=", "), "  [", nrow(x$extra), "x", ncol(x$df), "]\n", sep="")
  cat("Clusters: ", paste(names(x$clusters), collapse=", "), "\n", sep="")
}


#' Convert clusterfly object to data.frame.
#' Concatenates data and cluster assignments into one data.frame.
#' Cluster assignments are prefixed with \code{cl_}.
#'
#' @export
#' @method as.data.frame clusterfly
#' @param x clusterfly object
#' @param ... ignored
#' @keywords manip
as.data.frame.clusterfly <- function(x, ...) {
  cl <- as.data.frame(x$clusters)
  if (ncol(cl) > 0) {
    names(cl) <- paste("cl_", names(cl), sep="")
  } else {
    cl <- NULL
  }
  do.call("cbind", compact(list(x$df, x$extra, cl)))
}

#' Match all cluster indices to common reference.
#'
#' It's a good idea to run this before running any
#' animation sequences so that unnecessary colour
#' changes are minimised.
#'
#' @param cf clusterfly object
#' @param reference index to reference clustering
#' @param method method to use, see \code{\link{clarify}}
#' @keywords manip
#' @export
#' @examples
#' o <- olive_example()
#' o <- cfly_clarify(o, "Region")
cfly_clarify <- function(cf, reference=1, method="rowmax") {
  ref <- cf$clusters[[reference]]
  cf$clusters <- sapply(cf$cluster, function(x) clarify(x, ref, method=method), simplify=FALSE)
  cf
}

#' Add clustering.
#'
#' Clustering method needs to respond to \code{\link{clusters}},
#' if the default does not work, you will need to write
#' your own to extract clusters.
#'
#' @param cf clusterfly object
#' @param method clusterfing method (function)
#' @param ... arguments passed to clustering method
#' @param name name of clustering
#' @keywords manip
#' @export
#' @examples
#' o <- olive_example()
#' cfly_cluster(o, kmeans, 4)
#' cfly_cluster(o, kmeans, 4, name="blah")
cfly_cluster <- function(cf, method, ..., name = deparse(substitute(method))) {
  cf[[name]] <- clusters(method(cf$df, ...))
  cf
}


#' Dynamic plot: Animate glyph colours
#'
#' This function will animate until you manually break the loop
#' using Ctrl-Break or Ctrl-C.
#'
#' @param cf list of cluster ids that you want to animate across
#' @param clusters clusters to display
#' @param pause clusters number of seconds to pause between each change
#' @param print print current cluster to screen?
#' @param max_iterations maximum number of interations
#' @keywords dynamic
#' @export
#' @examples
#' # Press Ctrl-Break or Ctrl-C to exit
#' if (interactive()) {
#' o <- olive_example()
#' cfly_animate(cfly_clarify(o), max = 5)
#' close(o)
#' }
cfly_animate <- function(cf, clusters = seq_along(cf$clusters), pause = 1, print=TRUE, max_iterations = 100) {
  g <- cf$ggobi()
  gd <- g[1]

  if (is.character(clusters)) clusters <- match(clusters, names(cf$clusters))

  count <- 1
  while(TRUE) {
    for(i in clusters) {
      if (!valid_ggobi(g)) return()
      if (print) cat("Current cluster: ", names(cf$clusters)[i], "\n")
      glyph_colour(gd) <- cf$clusters[[i]]
      Sys.sleep(pause)

      count <- count + 1
      if (count > max_iterations) return()
    }
  }
}

