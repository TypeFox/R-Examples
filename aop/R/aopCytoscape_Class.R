#' aop_cytoscape class
#' 
#' Creates an object of class aop_cytoscape
#' 
#'  @section Slots:
#'    \describe{
#'      \item{\code{name}:}{Object of class \code{"character"}, containing the
#'      name of the AOP.}
#'      \item{\code{nodes}:}{Object of class \code{"list"}, containing the list
#'      of nodes.}
#'      \item{\code{edges}:}{Object of class \code{"list"}, contianing the list of
#'      edges.}
#'    }
#'  @name aop_cytoscape-class
#'  @rdname aop_cytoscape-class
#'  @aliases aop_cytoscape
#'  @exportClass aop_cytoscape
#'  @author Lyle D. Burgoon
#'  
#'  
#'  @export
aop_cytoscape <- setClass(
  #Set the name
  "aop_cytoscape",
  
  #Define the slots
  slots = c(
    name = "character",
    nodes = "list",
    edges = "list"
  ),
  
  validity=function(object){
    if(!is.character(object@name)){
      return("name must be character")
    }
    return(TRUE)
  }
  
)