#' Converts population graph to html file for interactive viewing
#' 
#' This function takes a population graph that has been 'decorated' with
#'  sufficient spatial data to make a html file that uses the D3 visualization
#'  javascript mateirals to view it interactively.
#' @param graph A \code{popgraph} object.
#' @param file The path to the html file to be saved.  If not given then the html
#'  text is returned by the function.
#' @return The text of the html file to be saved or viewed in the appropriate browser.
#' @author Rodney J. Dyer <rjdyer@@vcu.edu>
#' @export
to_html <- function( graph, file ) {
  if( !inherits( graph, "popgraph") )
    stop("Cannot save a html file from a popgraph that is not made from a popgraph...")
  
  heading <- system.file("extdata","d3header.html",package="popgraph")
  footing <- system.file("extdata","d3footer.html",package="popgraph")
  
  if( !nchar(heading) | !nchar(footing) )
    stop("Cannot run this until the package is actually installed. ")  
  
  head <- paste( readLines(heading),collapse="\n" )
  foot <- paste( readLines(footing), collapse="\n" )
  json <- to_json( graph )
  
  htmltext <- paste(head,json,foot,collapse="\n")
  
  if( !missing(file) ) {
    write(htmltext,file)
    invisible(htmltext)
  }
  
  else{
    return( htmltext)
  }
}