#' Convience function for file exports
#' 
#' This function is a chokepoint for exporting
#'  \code{popgraph} objects to other formats.
#' @param graph An object of type \code{popgraph}.
#' @param file The path to save the graph into.
#' @param format The type of output file to use.  Options are:
#'  \itemize{
#'    \item{json } { Save as JSON format}
#'    \item{kml } { Save to KML format to view in GoogleEarth}
#'    \item{graphml } { Save as Graph Markup Language}
#'    \item{html } { Save to an interactive html format viewable in your browser}
#'    \item{pajek } { Save as input to Pajek}
#'    \item{pgraph } { Save as input for GeneticStudio (default)}
#'    \item{adjacency } { Saves as an adjacency matrix in csv format}
#'    \item{paths } { Saves as shortest paths matrix in csv format}
#'    \item{weights } { Saves as weighted adjacency matrix in csv format}
#'  }
#' @param ... Ignored
#' @return Nothing 
#' @export
#' @author Rodney J. Dyer <rjdyer@@vcu.edu>
write_popgraph <- function(graph,file,format="pgraph",...){
  
  if(!is(graph,"popgraph"))
    stop("This requires a popgraph object")
  
  if( !(format %in% c("json","kml","graphml","html","pajek","pgraph","adjacency","paths","weights")))
    stop("Unrecognized output format.")
  
  if( missing(file) )
    stop("You need to pass a file to this function.")
  
  switch( format, 
          json=to_json(graph,file),
          html=to_html(graph,file),
          kml=to_kml(graph,file),
          pgraph=to_pgraph(graph,file),
          adjacency=function(graph,file){ 
            a <- to_matrix(graph,mode="adjacency")
            write.csv(a,file=file)
          },
          paths=function(graph,file){
            a <- to_matrix(graph,mode="shortest path")
            write.csv(a,file=file)
          },
          weights=function(graph,file){
            a <- to_matrix(graph,mode="edge weights")
          },
          write.graph(graph,file,format=format)
          )  
}