#' Riverplot styles
#'
#' Riverplot styles
#'
#' Riverplot styles are just lists with key-value pairs that define how
#' nodes and edges are drawn. Although there are attributes that are only
#' applicable to either nodes or edges, there are no separate style lists for
#' these objects. 
#' 
#' The \code{default.style} function simply returns the default style
#' defined in the riverplot package (including edge and node attributes).
#' 
#' The \code{updateRiverplotStyle} function updates all missing fields in
#' the \code{style} object with the styles from the \code{master} style.
#' 
#' When a node is drawn, the styles are determined by precedence. Command
#' line arguments to \code{\link{riverplot}()} function override any defined
#' styles. For all other parameters styles associated with nodes are used, and
#' if absent, inserted from the \code{default.style} argument to the
#' \code{\link{riverplot}()} function. If this argument is missing, style is
#' taken from the argument returned by the \code{default.style} function.
#'
#' Not recognized fields and values will be silently ignored.
#'
#' Following style fields and values are defined:
#'
#' \describe{
#' \item{nodestyle}{(default: regular). Values:
#'    \describe{
#'          \item{regular}{rectangular box with a label}
#'          \item{point}{a color dot}
#'          \item{invisible}{No node is drawn. This is used to seamlessly
#'                           integrate edges.}
#'    }}
#' \item{edgestyle}{(default: sin). Describes how the edge looks like.
#'    \describe{
#'         \item{sin}{A sinusoidal edge}
#'         \item{straight}{A straight edge}
#'    }
#'}
#' \item{edgecol}{(default: "gradient"). How edge color is generated.  Values:
#'    \describe{
#'         \item{gradient}{A color gradient generated based on parent and child
#' node that form the edge}
#'         \item{col}{The color specified in the "col" attribute of the edge}
#'    }}
#' \item{col}{(default: "grey"). Color of the node or edge (for edges, it
#' is used only if the "edgecol" attribute is "col".}
#' \item{srt}{(default: "90"). Rotation of the label (see \code{\link{par}})}
#' \item{lty}{(default: 1). Line type to draw around node and edges}
#' \item{textcol}{(default: "black"). Color of the node label.}
#' }
#' @param style style to update
#' @param master master style to use for updating
#' @return Both functions return an object of the riverplotStyle class
#' (which is, in fact, just a list with key-value pairs that you can access, inspect and manipulate
#' manually at will).
#' @author January Weiner
#' @examples
#' # To view the default style specification, type
#' default.style()
#' 
#' ex <- riverplot.example()
#' ds <- default.style()
#' plot( ex, default_style= ds )
#'
#' # nodes with unspecified style will now be semi-transparent red:
#' ds[["col"]] <- "#FF000099"
#' plot( ex, default_style= ds )
#' 
#' @name riverplot-styles
NULL
