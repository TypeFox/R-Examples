data.loadTemperature <- function(year, temperature.year, 
				temperature.previous.year, 
				from.previous.year.doy, length, 
				position, scale.factor=0.1){

	previous.year.length <- ifelse(util.isLeapYear(year-1),366,365)
	temperature.vec <- vector(mode="numeric", length=length)
	
	previous.year.count <- length(from.previous.year.doy:previous.year.length)
	temperature.vec[1:previous.year.count] <- 
		temperature.previous.year[from.previous.year.doy:previous.year.length,position]
	temperature.vec[(previous.year.count+1):length(temperature.vec)] <- 
		temperature.year[1:length((previous.year.count+1):length(temperature.vec)),position]

	return(temperature.vec * scale.factor)
}