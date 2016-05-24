# Converts Julian day of year (day of year, year)   
# into a date (day,month,year). If year is missing a non-leap year
# is assumed.
jul2date2 <- function(d,y) {
	if(missing(y)) y <- 2001
	
	res <- .C("Cjul2date2",doy=as.integer(d),year=as.integer(y),day=integer(1),month=integer(1),PACKAGE="pheno")

	return(list(day=res$day,month=res$month,year=res$year))
}
