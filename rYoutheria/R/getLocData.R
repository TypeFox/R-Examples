#' Get location information from YouTheria
#' 
#' Retrieves location information stored on YouTheria
#' 
#' @param country Character specifying the country within which to search for locations.
#' for a list of countries used getCountries().
#' @param StudyUnitId Numeric specifying the StudyUnitId to search for
#' @return A dataframe in which each rows gives the details of a study unit
#' @export        
#' @examples
#' \dontrun{
#' # Get a dataframe of Indian study units
#' Indian_StudyUnits <- getLocData(country = 'India')
#' }
#' 
getLocData <-
  function(country=NULL,StudyUnitId=NULL){
       
    if(!is.null(StudyUnitId)&!is.null(country)) stop('Cannot use both StudyUnitId and country at the same time')
    
    out <- runURL(paste('?id=', StudyUnitId,
                 '&country=', country, sep=''), 'l')
    
    if(length(out)==0){
      if(!is.null(StudyUnitId)) warning('No Data returned for this StudyUnitId(s).')
      if(!is.null(country)) warning('No Data returned for this country(s). Ensure you have capitalised appropriately')
    }
   
    return(out)
  }
