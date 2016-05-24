NULL
#' @name meteofranceSynop
#' @rdname meteofranceSynop
#' @title SYNOP Weather Data in France
#' @description \code{data.frame} table containing SYNOP weather data obtained through \code{MeteoDataFrance} R package (\url{https://github.com/ecor/MeteoDataFrance})
#' @docType data
#' @usage data(meteofranceSynop)
#' @format data.frame
#' @source \url{https://donneespubliques.meteofrance.fr}(in pariticular \url{https://donneespubliques.meteofrance.fr/?fond=produit&id_produit=94&id_rubrique=32})
#' @author Emanuele Cordano
#' 
#' 
#' @examples  
#' 
#' library(ggmap)
#' data(meteofranceSynop)
#' 
#' 
#' 
#' 
#' data=meteofranceSynop[meteofranceSynop$timestamp==meteofranceSynop$timestamp[250],]
#' 
#' 
#' 
#' map <- get_map(location ="France", zoom = 6)
#' 
#' size <- 3
#' 
#' gsynop <- ggmap(map) +
#' 		geom_point(data = data,aes(x = longitude, y = latitude),size=size,  alpha
#' 						=1, color="blue",show.legend  = FALSE)
#' 
#' ## Uncomment if you want to save in PDF format the otput of gsynop
#' ## ggsave("test-map.pdf", gsynop,width=10,height=10)
#'  
NULL