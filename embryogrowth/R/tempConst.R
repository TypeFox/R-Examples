#' Timeseries of constant temperatures for nests
#' @title Timeseries of constant temperatures for nests
#' @author Marc Girondot \email{marc.girondot@@u-psud.fr}
#' @docType data
#' @name tempConst
#' @description Timeseries of temperatures for nests
#' @references Girondot, M. & Kaska, Y. Submitted. A model to predict 
#'             temperature dependency on embryo growth rate and incubation
#'             duration from field data.
#' @keywords datasets
#' @usage tempConst
#' @examples
#' \dontrun{
#' library(embryogrowth)
#' # Same as:
#' # GenerateConstInc(durations = rep(104*60*24, 11),
#' # temperatures = 25:35,
#' # names = paste0("T",25:35))
#' data(tempConst)
#' tempConst_f <- FormatNests(tempConst)
#' x <- structure(c(118.768297442004, 475.750095909406, 306.243694918151, 
#' 116.055824800264), .Names = c("DHA", "DHH", "T12H", "Rho25"))
#' # pfixed <- c(K=82.33) or rK=82.33/39.33
#' pfixed <- c(rK=2.093313)
#' resultNest_4p <- searchR(parameters=x, fixed.parameters=pfixed, 
#' 	temperatures=formated, derivate=dydt.linear, M0=1.7, 
#' 	test=c(Mean=39.33, SD=1.92))
#' data(resultNest_4p)
#' # Use the fited parameters from resultNest_4p with  
#' # the constant incubation temperatures:
#' plot(resultNest_4p, temperatures=tempConst_f,  
#' 	stopattest=TRUE, series="all", xlim=c(0,120),  
#' 	ylimT=c(22, 32), show.stages=FALSE, show.PT=FALSE,  
#' 	show.temperatures=FALSE, show.TSP=FALSE)
#' }
#' @format A dataframe with raw data.
NULL
