yearlength <- function(year){
	year_length <- vector(mode="numeric", length=length(year) )

	leaps <- leapYears(year)
	year_length[leaps] <- 366
	year_length[!leaps] <- 365

	return(year_length)
}
