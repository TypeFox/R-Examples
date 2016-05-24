#'@title get climate data names in GLTC dataset
#'@description get climate names for the Global Lake Temperature Collaboration dataset.
#'If a \code{lake_name} is used, only names of climate drivers that exist for that lake will be returned. 
#'If no \code{lake_name} is specified, all climate driver names for the entire dataset will be returned. 
#'
#'@param lake_name a valid name of a lake in the GLTC dataset (see \code{\link{get_lake_names}}).
#'@return a character vector of valid climate variable names
#'@importFrom dplyr filter
#'@seealso \code{\link{get_climate}}, \code{\link{get_lake_names}}
#'@examples
#'get_climate_names()
#'get_climate_names('Victoria')
#'@export
get_climate_names <- function(lake_name){
  skip_names <- c(temp_types() , 'X')
  
  # -- fix for R CMD check 'no visible binding for global variable'
  siteID <- "_private"
  # -- fix for R CMD check 'no visible binding for global variable'
  
  if (missing(lake_name)){
    
    val_names <- unique(gltc_values$variable)
    
  } else {
    data <- filter(gltc_values, siteID == get_site_ID(lake_name))
    val_names <- unique(data$variable)
  }
  
  climate_vars <- val_names[!val_names %in% skip_names]
  
  return(climate_vars)
}