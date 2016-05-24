#' converting epiJSON into a list of lists
#' 
#' takes an epiJSON string or file and converts to a list of lists  
#' later this needs to be made more formal


#' @param file an epiJSON filename or string to convert to R
#' 
#' @return a list of lists of the epijson content
#' @examples
#' listJSON <- epiJSON2r( system.file("extdata//example.JSON", package="repijson"))
#' str(listJSON)
#' #from within the package would do this
#' #listJSON <- epiJSON2r("extdata//example.JSON")
#'  
#' @export
#' 
epiJSON2r <- function(file)
{
  
  #first test output the string
  cat(file)
  
  ## convert from json to list
  #listJSON <- RJSONIO::fromJSON(file)
  listJSON <- jsonlite::fromJSON(file) 
  
  
  return(listJSON)
}
