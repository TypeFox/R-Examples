#' Plot an Aster graph object comprised of the vertice and edge tables.
#' 
#' Function for obtaining and plotting graph (network) objects using ggplot2 and ggnet2 functions.
#' 
#' @param channel connection object as returned by \code{\link{odbcConnect}}
#' @param graph an object of class \code{'toagraph'} referencing graph 
#'   tables in Aster database.
#' @param v a SQL \code{SELECT} that returns key values or a list of key values (corresponding 
#'   to the \code{vertex.names} attribute) of the vertices to include in the graph. 
#'   When not \code{NULL} this guarentees that no other vertices or edges between other vertices 
#'   are included in the resulting network.
#' @param vertexWhere optionally, a \code{SQL WHERE} clause to subset vertex table. When not \code{NULL}
#'   it overrides \code{vertexWhere} condition from the \code{graph}.
#' @param edgeWhere optionally, a \code{SQL WHERE} clause to subset edge table. When not \code{NULL}
#'   it overrides \code{edgeWhere} condition from the \code{graph}.
#' @param ... Other arguments passed on to \code{\link[GGally]{ggnet2}} visualization function. 
#' @param allTables pre-built information about existing tables.
#' @param test logical: if TRUE show what would be done, only (similar to parameter \code{test} in \pkg{RODBC} 
#'   functions: \link{sqlQuery} and \link{sqlSave}).
#'    
#' @return a ggplot object
#' @export
#' @examples 
#' if(interactive()) {
#' # initialize connection to Lahman baseball database in Aster 
#' conn = odbcDriverConnect(connection="driver={Aster ODBC Driver};
#'                          server=<dbhost>;port=2406;database=<dbname>;uid=<user>;pwd=<pw>")
#' 
#' policeGraphUn = toaGraph("dallaspolice_officer_vertices", "dallaspolice_officer_edges_un", 
#'                          directed = FALSE, key = "officer", 
#'                          source = "officer1", target = "officer2", 
#'                          vertexAttrnames = c("offense_count"), edgeAttrnames = c("weight"))
#' 
#' # visualize the whole graph                         
#' showGraph(conn, policeGraphUn, node.label = "vertex.names", node.size="offense_count")
#' 
#' # visualize subgraph using both vertex and edge filters
#' showGraph(conn, policeGraphUn, 
#'           vertexWhere = "officer ~ '[A-Z ].*'", edgeWhere = "weight > 0.10", 
#'           node.label = "vertex.names", node.size="offense_count")
#' 
#' }
showGraph <- function(channel, graph, v=NULL, 
                      vertexWhere=graph$vertexWhere, edgeWhere=graph$edgeWhere,
                      ..., allTables=NULL, test=FALSE) {
  
  if (test && is.null(allTables))
    stop("Must provide allTables when test==TRUE.")
  
  net = computeGraph(channel, graph, v, vertexWhere = vertexWhere, edgeWhere = edgeWhere,
                     allTables = allTables, test = test)
  
  if (test)
    return(net)
  
  p = GGally::ggnet2(net, ...) 
  
  return(p)
}