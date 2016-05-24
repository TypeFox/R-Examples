# Summarize water balance fluxes, compare to simulated change in soil water storage.
#' Checks MAESPA water balance by adding up fluxes for the soil water balance calculated by MAESPA. Prints the
#' total fluxes (precipitation, transpiration, canopy interception, and so on),
#' and the difference between ingoing and outgoing. This is a
#' debugging/checking tool: if the model has serious problems, there will be a
#' missing sink or source.
#' 
#' @param x A dataframe returned by \code{\link{readwatbal}}.
#' @param usemeaset Whether to use measured ET from the met file, or simulated
#' (Default).
#' @author Remko Duursma
#' @seealso \code{\link{readwatbal}}
#' @keywords misc
#' @export
checkwatbal <- function(x = readwatbal(), usemeaset=FALSE){
	
	if(is.na(x$wsoilroot[1]))
		startat <- 2
	else
		startat <- 1
	
	# delta storage
	n <- nrow(x)
	delstor <- x$wsoil[n] - x$wsoil[1]

	# components of water balance.
	tfall <- sum(x$tfall[startat:n])
	etmm <- sum(x$et[startat:n])
	overflow <- sum(x$overflow[startat:n])
	
	if(usemeaset){
		et <- x$etmeas[startat:n]
		et[et < 0] <- 0
		etmm <- sum(et)
	}

	discharge <- sum(x$discharge[startat:n])
	soilevap <- sum(x$soilevap[startat:n])
	ppt <- sum(x$ppt[startat:n])
	
	evapstore <- sum(x$evapstore[startat:n])
	
	# sum of fluxes (= delta storage??)
	delflux <- tfall - etmm - discharge - soilevap - overflow

	cat("Precip:      ", ppt, "\n")
	cat("Throughfall: ", tfall, "\n")
	cat("Wet ET:      ", evapstore, "\n\n")
	cat("ET:          ", etmm, "\n")
	cat("Discharge:   ", discharge, "\n")
	cat("Overflow:    ", overflow, "\n")
	cat("Soil evap.:  ",soilevap, "\n")
	cat("\n")
	cat("Change in storage:", delstor, "\n")
	cat("Total fluxes:     ", delflux, "\n")
	cat("\n")
	cat("Error of", delstor-delflux, "mm \n")
	
}
