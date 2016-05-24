#' @name minard
#' @title Minard Napoleon Russian campaign data
#' @description The data set used by Charles Joseph Minard to generate the
#' famous graph. The example below shows how to recreate the
#' main panel of the graph using riverplot from the provided data.
#' 
#' First, node and edge data frames must get new column names (see
#' \code{\link{makeRiver}} function for details). Then, based on the
#' direction of the Napoleon army, style information (right and left edge
#' color style for each node) is entered in the
#' \var{nodes} variable. Then, a riverplot object is generated from the
#' nodes and edges data frames.
#'
#' To use the same color coding as Minard, the \var{direction} variable is
#' converted to color codes in the \var{col} column of the
#' \var{edges} object.
#'
#' Finally, a plot is created using \code{lty=1} and a style in which nodes
#' are not shown, and the edges are straight (like in the original Minard
#' plot) rather than curved.
#'
#' @docType data
#' @usage minard
#' @format Named list with two data frames: 
#' \describe{
#' \item{nodes}{data frame with geographic locations of the Napoleon army
#' (longitude and latitude) and the direction of the march}
#' \item{edges}{connections between positions}
#'}
#' @source Charles Joseph Minard
#' @examples
#' data( minard )
#' nodes <- minard$nodes
#' edges <- minard$edges
#' colnames( nodes ) <- c( "ID", "x", "y" )
#' colnames( edges ) <- c( "N1", "N2", "Value", "direction" )
#'
#' # color the edges by troop movement direction
#' edges$col <- c( "#e5cbaa", "black" )[ factor( edges$direction ) ]
#' # color edges by their color rather than by gradient between the nodes
#' edges$edgecol <- "col"
#'
#' # generate the riverplot object and a style
#' river <- makeRiver( nodes, edges )
#' style <- list( edgestyle= "straight", nodestyle= "invisible" )
#'
#' # plot the generated object
#' plot( river, lty= 1, default_style= style )
#' # Add cities
#' with( minard$cities, points( Longitude, Latitude, pch= 19 ) )
#' with( minard$cities, text( Longitude, Latitude, Name, adj= c( 0, 0 ) ) )
#'
#' @author January Weiner
NULL
