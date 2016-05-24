# Converts Julian day of year (day of year, year)   
# into a string date "DD.MM.YYYY". If year is missing a non-leap year
# is assumed.
jul2date1 <- function(d,y) {
	if(missing(y)) y <- 2001
	
	res <- .C("Cjul2date1",doy=as.integer(d),year=as.integer(y),date="DD.MM.YYYY",PACKAGE="pheno")

	return(res$date)
}
