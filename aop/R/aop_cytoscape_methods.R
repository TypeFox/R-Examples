#' @include aopCytoscape_Class.R
NULL

#' Get Node Name from ID
#' 
#' 
#' Given an id, this method returns an \code{aop_cytoscape} node name.
#' 
#' @param theObject is an AOP as an object of class \code{aop_cytoscape}.
#' 
#' @param id an object of class \code{character} such as "389".
#' 
#' @return the name of the node
#' 
#' @export
#' @docType methods
#' @rdname aop_cytoscape-methods
#' 
#' @examples
#' library(graph)
#' steatosis_json_file <- system.file("extdata", "steatosis_aop_json.cyjs",
#' package = "aop")
#' steatosis_aop <- convert_cytoscape_to_aop(steatosis_json_file)
#' getAOPNodeName(steatosis_aop, "389") 
setGeneric(name="getAOPNodeName", 
           def=function(theObject, id){
             standardGeneric("getAOPNodeName")
           }
           )

#' @rdname aop_cytoscape-methods
setMethod(f="getAOPNodeName",
          signature="aop_cytoscape",
          definition=function(theObject, id){
            return(theObject@nodes[[id]]$data$name)
          }
          )