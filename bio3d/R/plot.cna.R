plot.cna <- function(x, pdb=NULL, weights=NULL, vertex.size=NULL,
                     layout=NULL, col=NULL, full=FALSE, scale = TRUE, color.edge = FALSE, ...) {

  ##- Function for plotting cna networks the way we like them.
  ##   Returns the plot layout coordinates silently. These can 
  ##   be passed to identify.cna() 
  ##
  ##- Examples:
  ##   plot.cna(net)
  ##   plot.cna(net, pdb)
  ##   plot.cna(net, layout=layout.cna(net,pdb))
  ##   plot.cna(net, layout=layout.cna(net,pdb, full=TRUE), full=TRUE)
  ##   plot.cna(net, full=T, layout=layout.cna(net,pdb, full=T), vertex.size=3, weights=1, vertex.label=NA)
  ##
  ##- Other options:
  ##   \dots can contain all ?igraph.plotting options, including:
  ##   col=vmd.colors(),
  ##   mark.groups=list() - A list of vertex id vectors 
  ##   mark.col=vmd.colors(alpha=0.3),
  ##   mark.border=vmd.colors()
  ##     etc... see ?plot.igraph
  ##  AND
  ##   vertex.size:  Node sizes:   V(x$network)$size
  ##   vertex.color: Node colors:  V(x$network)$color
  ##   vertex.label: Node labels:  V(x$network)$name - use NA to omit
  ##   edge.width:   Edge weights: E(x$network)$weight
  ##   edge.color:   Edge colors:  E(x$network)$color
  ##    (also vertex.label.color, vertex.label.cex etc.
  ##     see ?igraph.plotting)
  ##   

  ## Check for presence of igraph package
  oops <- requireNamespace("igraph", quietly = TRUE)
  if (!oops) {
     stop("igraph package missing: Please install, see: ?install.packages")
  }
  
#  if(color.edge) {
#     oops <- require(classInt)
#     if (!oops) {
#        warning("package classInt missing: color.edge is set to FALSE.
#            To make color.edge work, please install the missing package. See: ?install.packages")
#        color.edge = FALSE
#     }    
#  }

  ##- Determine which network to plot along with node size
  if(full) {
    ## Plot the 'full' all-atom network
    y <- x$network 

    if(is.null(vertex.size)) {
      ## Scale up the node size to something visible
      if(max(igraph::V(y)$size) < 10) {
        igraph::V(y)$size = igraph::V(y)$size + 13
      }
    }
  } else {
    ## Plot the 'coarse-grained' community network
    y <- x$community.network
    ## here we will leave the node size as is
  }

 
  ##- Determine edge weights and scale values for plotting
  if(is.null(weights)){
    ## Use weights defined in network
    weights <- igraph::E(y)$weight
    
    if(is.null(x$call$minus.log)){  
      ## If '$call$mins.log' is NULL => -log option was used in cna()
      ##  so we we will revert back with exponential here
      weights <- exp(-weights)
    } else{
      if(x$call$minus.log){
        ## Again here we have 'minus.log=TRUE'
        weights <- exp(-weights)
      }
    }
    ## Lets scale the weights to lie between 1 and 5
#    weights <- (weights - min(weights)) / max(weights - min(weights)) * (1 - 5) + 5
    if(scale && (length(weights)>1)) weights <- (weights - min(weights)) / max(weights - min(weights)) * 4 + 1
    else weights <- 10 * weights
  }
  
  ##- Obtain the plot layout coords
  if(!is.null(pdb) && is.null(layout)) {
    cat("Obtaining layout from PDB structure\n")
    layout = layout.cna(x, pdb, full=full)
  }
  if(is.null(pdb) && is.null(layout)) {
    cat("Obtaining estimated layout with fruchterman.reingold\n")
    layout <- igraph::layout.fruchterman.reingold(y, weights=weights)
  }
  if(dim(layout)[2] != 2){
    stop("Input 'layout' must be an Nx2 matrix, where N is the number of communities")
  }
  
  if(color.edge) {

#     vec2color <- function(vec, pal=c("blue", "green", "red"), n=10) {
#        ##-- Define a color scale from a numeric vector
#        return( findColours(classIntervals(vec, n=n, style="equal"), pal) )
#     }
     vec2color <- function(vec, pal=c("blue", "green", "red"), n=30) {
        col <- colorRampPalette(pal)(n)
        vec.cut <- cut(vec, seq(min(vec), max(vec), length.out=n), include.lowest = TRUE)
        levels(vec.cut) <- 1:length(col)
        return( col[vec.cut] )
     }

     colors <- vec2color(weights)
     igraph::plot.igraph(y, edge.width=weights, edge.color = colors, layout=layout, vertex.color=col, vertex.size=vertex.size, ...)

  } else { 

     igraph::plot.igraph(y, edge.width=weights, layout=layout, vertex.color=col, vertex.size=vertex.size, ...)

  }
  
  ## Silently return plot coordinates
  #class(layout) = "cna"
  layout <- layout
}
