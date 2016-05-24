NULL
#' @name meteofrance
#' @rdname meteofrance
#' @title Snow Weather Data in France
#' 
#' @docType data
#' @usage data(meteofrance)
#' @format data.frame
#' @source \url{https://donneespubliques.meteofrance.fr}(in pariticular \url{https://donneespubliques.meteofrance.fr/?fond=produit&id_produit=94&id_rubrique=32})
#' @author Emanuele Cordano
#' 
#' 
#' @description The \code{meteofrance} data frame is a  table containing snow and weather data obtained through \code{MeteoDataFrance} R package (\url{https://github.com/ecor/MeteoDataFrance}).
#' @details It is organized with the following fields: 
#' \describe{
#'  \item{\code{"latitude"}:}{ latitude of the weather station; } 
#'  \item{\code{"longitude"}:}{ longitude of the weather station; }
#'  \item{\code{"station_id"}:}{ id of the weather station;}         
#'  \item{\code{"altitude"}:}{  altitude of the weather station;}
#'  \item{\code{"location"}:}{location name  of the weather station;}         
#'  \item{\code{"timestamp"}:}{\code{\link{POSIXlt}} date and time:}             
#'  \item{\code{"VW","TA","TD","RH",...}:}{ weather variables(name in accordance with SMET specification);}             
#'  \item{\code{"*_reapeted"}:}{ the column is repeated};
#' }
#' The \code{meteofrance} data frem has an attribute: \code{metaparam} is a data frame containing meta info on the weather variables:
#' \describe{  
#' \item{\code{"Descriptif"}:}{ weather variable description provided by MeteoFrance (in French);}     
#' \item{\code{"IDparam","Mnemonique"}:}{ MeteoFrance id for the weather variable;}    
#' \item{\code{"type"}:}{ weather variable type;} 
#' \item{\code{"unite"}:}{ measerement unit of the weather variable;}           
#' \item{\code{"SMET_ID"}:}{ SMET id for the weather variable (they are used as field/column names in \code{meteofrance});}    
#' \item{\code{"SMET_UNIT_MULTIPIPLIER"}:}{  SMET unit multiplier respect to the SMET_ID variable's MKSA unit (see SMET specifications);}      
#' \item{\code{"SMET_UNIT_OFFSET"}:}{ SMET unit offset respect to the SMET_ID variable's MKSA unit (see SMET specifications).}      
#' 
#' }
#'  
#' 
#' 
#' 
#' 
#' 
#' @seealso \code{\link{as.smet}}

#' 
#' @examples  
#' 
#' library(ggmap)
#' data(meteofrance)
#' 
#' dates <- as.Date(meteofrance$timestamp)
#' 
#' data=meteofrance[dates==dates[1],]
#' 
#' 
#' 
#' map <- get_map(location ="France", zoom = 6)
#' 
#' size <- 3
#' 
#' gsnow <- ggmap(map) +
#' 		geom_point(data = data,aes(x = longitude, y = latitude),size=size,  alpha
#' 						=1, color="blue",show.legend  = FALSE)
#' 
#' ## Uncomment if you want to save in PDF format the otput of gsnow
#' ## ggsave("test-map.pdf", gsnow,width=10,height=10)
#'  
NULL