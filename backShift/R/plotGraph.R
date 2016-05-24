#' Plotting function to visualize directed graphs
#'
#' @description Given a point estimate of the connectivety matrix or the
#'  adjacency matrix, this function visualizes the directed graph using
#'  \code{\link[igraph]{plot.igraph}} from the package \code{\link[igraph]{igraph}}. If a point estimate
#'  is plotted, the edges' intensity reflects the magnitude of the coefficients.
#'  If the result is an adjacency matrix estimated by stability selection then
#'  the edges' width reflects how often an edge was selected and the intensity
#'  reflects the magnitude of the coefficients (if this information is also
#'  provided).
#'
#' @details Currently not all options of \code{\link[igraph]{igraph}} are used; additional
#'  arguments are ignored.
#'
#' @param estimate Estimate of connectivity matrix. This can be a point estimate
#'  with entry \eqn{A_{ij}} being the estimated edge weight for the edge from
#'  node \eqn{i} to node \eqn{j}. Otherwise, it can be the estimated adjacency
#'  matrix by a stability selection procedure as in \code{\link{backShift}}. In
#'  this case, the entry \eqn{A_{ij}} indicates how often the edge from
#'  node \eqn{i} to node \eqn{j} was selected.
#'  @param plotStabSelec Set to TRUE if \code{estimate} results from the stability
#'  selection procedure. Otherwise, \code{estimate} is assumed to be a point
#'  estimate.
#' @param labels Variable labels to be displayed in plot.
#' @param thres.point Value at which the point estimate should be thresholded,
#'  i.e. edges with coefficients smaller than \code{thres.point} are not
#'  displayed.
#' @param edgeWeights If stability selection result should be visualized,
#'  provide edgeWeights as a (pxp)-matrix to display the magnitude of the
#'  coefficients as the intensity of the edges.
#' @param thres.stab Indicate the threhold value that was used in the stability
#'  selection procedure. Used to determine the width of the plotted edges.
#' @param main Provide the title of the plot.
#' @param edge.color Color of the edges. Defaults to blue.
#' @param ... Optional arguments passed to the plotting function. Consists of
#'  igraph-type options like vertex.label.cex,vertex.label.color, edge.arrow.size
#'  or vertex.size etc.
#'
#'  @examples
#' # create a matrix A to be visualized
#' p <- 3
#' A <- diag(p)*0
#' A[1,2] <- 0.8
#' A[2,3] <- -0.8
#' A[3,1] <- 0.8
#'
#' # add column names to use as labels for nodes
#' colnames(A) <- c("1", "2", "3")
#'
#' # plot
#' plotGraphEdgeAttr(estimate = A, plotStabSelec = FALSE, labels = colnames(A),
#'                    thres.point = 0, thres.stab = NULL, main = "True graph")
plotGraphEdgeAttr <- function(estimate, plotStabSelec, labels, thres.point,
                              edgeWeights = NULL, thres.stab = 0.75,
                              main = "", edge.color = "blue", ...){

  # labels
  colnames(estimate) <- rownames(estimate) <- labels
  
  # estimate
  df.A <- melt(as.matrix(estimate))
  colnames(df.A) <- c("from", "to", "edge")

  is.empty <- sum(estimate) == 0
  if(plotStabSelec & !is.empty){
    # are edge weights provided when stability selection result
    # should be visualized?
    if(is.null(edgeWeights)){
      edgeWeights <- estimate
      edgeWeights[abs(edgeWeights) > 0] <- 1
    }
    df.ew <- melt(as.matrix(edgeWeights))
    colnames(df.ew) <- c("from", "to", "edge")

    combined <- df.A
    combined$estimate <- df.ew$edge
    combined.wo.zeros <- combined[-which(combined$edge == 0),]

    if(nrow(combined.wo.zeros) == 0){
        is.empty <- TRUE
    }else{
      # draw graph
      bio.network <- graph.data.frame(combined.wo.zeros,
                                      directed=TRUE,
                                      vertices = colnames(estimate))
      layoutfunction <-  layout.circle
      layout <- layoutfunction(bio.network)

      # using the absolute value of the point estimate for the
      # transparency of the chosen color
      color.attribute <- abs(E(bio.network)$estimate)
      alphas.to.use <- convert.given.min(color.attribute, 0,
                                         max(color.attribute),
                                         range.min = 0.2, range.max = 0.8)
      mycolors <- vary.alpha(alphas.to.use, edge.color)
      E(bio.network)$color <- mycolors

      # using the number of times selected in the stability selection as the width
      width.attribute <- abs(E(bio.network)$edge)
      E(bio.network)$width <- convert.given.min(width.attribute,
                                                thres.stab*100, 100,
                                                range.min = 3, range.max = 8)
    }

    }else if(!is.empty){
      combined <- df.A
      combined.wo.zeros <- combined[-which(abs(combined$edge) <= thres.point),]

      if(nrow(combined.wo.zeros) == 0){
        is.empty <- TRUE
        estimate[abs(estimate) <= thres.point] <- 0
      }else{
        # draw graph
        bio.network <-graph.data.frame(combined.wo.zeros,
                                       directed=TRUE,
                                       vertices = colnames(estimate))
        layoutfunction <- layout.circle
        layout <- layoutfunction(bio.network)

        # using the absolute value of the point estimate for
        # the transparency of the chosen color
        color.attribute <- abs(E(bio.network)$edge)
        rmin <- 0.2
        rmax <- 0.8
        alphas.to.use <- convert.given.min(color.attribute, 0,
                                           max(color.attribute),
                                           range.min = rmin, range.max = rmax)
        if(any(is.na(alphas.to.use)))
          alphas.to.use <- rep(rmax, length(color.attribute))
        mycolors <- vary.alpha(alphas.to.use, edge.color)
        E(bio.network)$color <- mycolors

        # using the number of times selected in the stability selection as the width
        width.attribute <- abs(E(bio.network)$edge)
        rmin.w <- 5
        rmax.w <- 5
        widths.to.use <- convert(width.attribute,
                                 range.min = rmin.w, range.max = rmax.w)
        if(any(is.na(widths.to.use)))
          widths.to.use <- rep(rmax.w, length(width.attribute))
        E(bio.network)$width <- widths.to.use
      }
    }

  optionals <- list(...)
  if(is.null(optionals$vertex.label.cex)) optionals$vertex.label.cex <-  2.5
  if(is.null(optionals$vertex.label.color)) optionals$vertex.label.color <-"black"
  if(is.null(optionals$vertex.color)) optionals$vertex.color <-"white"
  if(is.null(optionals$vertex.frame.color)) optionals$vertex.frame.color <-"black"
  if(is.null(optionals$edge.arrow.size)) optionals$edge.arrow.size <-1
  if(is.null(optionals$edge.arrow.width)) optionals$edge.arrow.width <- 1.5
  if(is.null(optionals$vertex.size)) optionals$vertex.size <-30
  if(is.null(optionals$vertex.label.dist)) optionals$vertex.label.dist <-0
  if(is.null(optionals$vertex.label.degree)) optionals$vertex.label.degree <- -pi/2

  if(!is.empty){
    plot(bio.network, layout=layout, main = main,
         vertex.shape="circle",
         vertex.label.cex=optionals$vertex.label.cex,
         vertex.label.color=optionals$vertex.label.color,
         vertex.color=optionals$vertex.color,
         vertex.frame.color=optionals$vertex.frame.color,
         vertex.size=optionals$vertex.size,
         vertex.label.dist=optionals$vertex.label.dist,
         vertex.label.degree=optionals$vertex.label.degree,
         edge.arrow.size = optionals$edge.arrow.size,
         edge.arrow.width = optionals$edge.arrow.width,
         edge.curved=autocurve.edges2(bio.network, start = 0.25))
  }else{
    plotGraph(estimate,main=main,labels=labels,
              vertex.shape="circle",
              vertex.label.cex=optionals$vertex.label.cex,
              vertex.label.color=optionals$vertex.label.color,
              vertex.color=optionals$vertex.color,
              vertex.frame.color=optionals$vertex.frame.color,
              vertex.size=optionals$vertex.size,
              vertex.label.dist=optionals$vertex.label.dist,
              vertex.label.degree=optionals$vertex.label.degree,
              edge.arrow.size = optionals$edge.arrow.size,
              edge.arrow.width = optionals$edge.arrow.width)
  }
}



plotGraph <- function(A, main="",labels=NULL,layoutfunction=layout.circle,...){
  if(is(A, "dgCMatrix")) A <- as.matrix(A)
  if(!is.matrix(A)) stop("A needs to be a matrix")
  if(nrow(A)!=ncol(A)) stop("A needs to have as many rows as columns")
  if(is.null(labels)) labels <-
      if( !is.null(cc <- colnames(A))) cc else as.character(1:ncol(A))
  G <- graph.adjacency(A,mode="directed",weighted="a")
  if(is.null(layoutfunction)) layoutfunction <-  layout.circle
  layout <- layoutfunction(G)

  optionals <- list(...)
  if(is.null(optionals$vertex.label.cex)) optionals$vertex.label.cex <- 1.5
  if(is.null(optionals$vertex.label.color))
    optionals$vertex.label.color <-rgb(0.8,0.1,0.1,0.7)
  if(is.null(optionals$vertex.color)) optionals$vertex.color <-"white"
  if(is.null(optionals$vertex.frame.color))
    optionals$vertex.frame.color <-rgb(0.8,0.1,0.1,0.5)
  if(is.null(optionals$edge.color)) optionals$edge.color <-rgb(0.1,0.1,0.1,0.5)
  if(is.null(optionals$edge.arrow.size)) optionals$edge.arrow.size <-0.7
  if(is.null(optionals$edge.arrow.width)) optionals$edge.arrow.width <-2
  if(is.null(optionals$vertex.size)) optionals$vertex.size <-30
  if(is.null(optionals$vertex.label.dist)) optionals$vertex.label.dist <-0
  if(is.null(optionals$vertex.label.degree)) optionals$vertex.label.degree <- -pi/2

  plot(G, layout=layout, main=main,
       vertex.label=labels,
       vertex.shape="circle",
       vertex.label.cex=optionals$vertex.label.cex,
       vertex.label.color=optionals$vertex.label.color,
       vertex.color=optionals$vertex.color,
       vertex.frame.color=optionals$vertex.frame.color,
       edge.color=optionals$edge.color,
       edge.arrow.size=optionals$edge.arrow.size,
       edge.arrow.width=optionals$edge.arrow.width,
       vertex.size=optionals$vertex.size,
       vertex.label.dist=optionals$vertex.label.dist,
       vertex.label.degree=optionals$vertex.label.degree)
}
