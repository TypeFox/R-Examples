#' Get Measurement Types
#' 
#' Retrieves a \code{data.frame} of measurement types available from YouTheria.
#' 
#' @param measurementType If \code{NULL} (default), then all measurement types are
#'        returned. Can also be 'numeric' or 'character' (or a list of either type)
#'        and will filter by Id and Name respectivly in the resulting data.frame.
#' @return A dataframe of measurement types giving their Id and Name        
#' @export
#' @examples
#' \dontrun{
#' # Get a dataframe of all measurement types
#' AllMT <- getMeasurementTypes()
#' 
#' # Seach by name
#' BM_MT <- getMeasurementTypes('Body Mass')
#' BM_LL_MT <- getMeasurementTypes(c('Body Mass','Limb Length'))
#' 
#' # Search by ID
#' MT1 <- getMeasurementTypes(1)
#' MT123 <- getMeasurementTypes(1:3)
#' }

getMeasurementTypes <-
function(measurementType = NULL){

  if(is.null(measurementType)){
    
    out <- runURL(type = 't')
    
  } else if(class(measurementType)=='character') {
    
    url<-paste('?MeasurementType=',measurementType,sep='')
    out <- runURL(URL = url, type = 't')
    
  } else if(class(measurementType) == 'numeric' | class(measurementType) == 'integer'){
    
    url<-paste('?id=',measurementType,sep='')
    out <- runURL(URL = url, type = 't')
    
  } else{
    
    stop('argument must be numeric, integer, charater or NULL')
    
  }
  
  if(is.null(out)){
    
    warning('no matching measurement types found')
    
  } else {
    
    out$Name <- as.character(out$Name)
    
  }
  
  return(out)
  
}
