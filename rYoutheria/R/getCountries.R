#' Get country names from YouTheria
#' 
#' Retrieves a \code{data.frame} of country names and IDs from YouTheria.
#' 
#' @return A dataframe of country names and IDs. These names can be used in 
#' \link{getMeasurementData} to restrict the search to a specific country
#' @export        
#' @examples
#' \dontrun{
#' # Get a dataframe of all countries
#' getCountries()
#' }

getCountries <- function(){
  
  listCountries <- runURL(type = 'c')
  return(listCountries)
  
}