#'Assign GBIF style degree cell ids
#'
#'Calculate and assign a GBIF-style degree cell id and centi-degree (0.1 
#'degrees, dividing a 1 degree cell into 100 centi-degree cells) cell id to each
#'record. This is a necessary previous step for some functions like
#'\code{\link{mapgrid}}
#'
#'@param indf input data frame containing biodiversity data set
#'@export
#'@family Data preparation functions
#'@examples \dontrun{
#'getcellid(inat)
#'}
getcellid <- function (indf){
  indf$Latitude <- as.numeric(indf$Latitude)
  indf$Longitude <- as.numeric(indf$Longitude)
  indf$Cell_id <- (((indf$Latitude %/% 1) + 90) * 360) + ((indf$Longitude %/% 1) + 180)
  indf$Centi_cell_id <- ((((indf$Latitude %% 1) * 10) %/% 1 ) * 10) + (((indf$Longitude %% 1) * 10) %/% 1)
  return(indf)
}
