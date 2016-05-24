
#' Benchmark experiment graph
#'
#' The benchmark summary plot takes the individual benchmark
#' experiment results into account. The y-axis represents the
#' data sets, the x-axis a podium with as many places as
#' candidate algorithms.
#'
#' @export
bsgraph0 <- function(x, ...) {
  UseMethod('bsgraph0')
}


#' @param x A \code{\link{dist}} or \code{\link[graph]{graphNEL-class}} object
#' @param ndists.show The number of distance levels to show
#' @param edge.col The color of edges (one or one for each distance level)
#' @param edge.lwd The line width of edges (one or one for each distance level)
#' @param node.fill The colors of nodes
#' @param ... Arguments passed to underlying function
#' @return The return value of \code{\link{bsgraph0.graphNEL}}
#' @method bsgraph0 dist
#' @family algperf-visualization
#' @S3method bsgraph0 dist
#' @rdname bsgraph0
bsgraph0.dist <- function(x, ndists.show = length(sort(unique(x))),
                          edge.col = gray(0.7), edge.lwd = 1,
                          node.fill = NULL, ...) {

  data <- as.matrix(x)

  nodes <- colnames(data)
  nnodes <- length(nodes)

  dists <- sort(unique(x))
  ndists <- length(dists)
  dshow <- dists[seq_len(ndists.show)]
  ndshow <- length(dshow)

  edge.col <- rep(edge.col, ndshow)
  edge.lwd <- rep(edge.lwd, ndshow)
  edge.len <- ceiling((1.2)^(seq_len(ndists)-1))
  edge.weight <- rev(seq_len(ndists))
  edge.lty <- c(rep('solid', ndshow),
                rep('blank', length(dists)-ndshow))

  graph <- new('graphNEL', nodes=nodes, edgemode='undirected')
  edgeAttrs <- list()
  nodeAttrs <- list()

  for ( i in 1:(nnodes-1) ) {
    for ( j in (i+1):nnodes ) {
      s <- data[i,j]

      if ( s %in% dshow ) {
        t <- which(s == dists)

        graph <- addEdge(nodes[i], nodes[j], graph, edge.weight[t])

        n <- paste(nodes[i], nodes[j], sep='~')
        edgeAttrs$len[n] <- edge.len[t]
        edgeAttrs$color[n] <- edge.col[t]
        edgeAttrs$lwd[n] <- edge.lwd[t]
        edgeAttrs$lty[n] <- edge.lty[t]
      }
    }
  }

  if ( !is.null(node.fill) )
    nodeAttrs$fillcolor[nodes] <- node.fill

  bsgraph0(graph, nodeAttrs=nodeAttrs, edgeAttrs=edgeAttrs)
}



#' @param layoutType Defines the layout engine
#' @method bsgraph0 graphNEL
#' @S3method bsgraph0 graphNEL
#' @rdname bsgraph0
bsgraph0.graphNEL <- function(x, layoutType = 'neato', ...) {

  attrs <- getDefaultAttrs(layoutType=layoutType)
  attrs$node$fixedsize <- TRUE
  attrs$node$fontsize <- 20

  ag <- agopen(x, '', layoutType=layoutType, attrs = attrs, ...)
  plot(ag)

  # Redraw nodes for beauty:
  par(new=TRUE)
  ag2 <- ag
  ag2@AgEdge <- list()
  plot(ag2)


  invisible(ag)
}


