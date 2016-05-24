# Converts a date (day,month,year) into a Julian date (doy,year)
date2jul2 <- function(d,m,y) {
	if(missing(y)) y <- 2000

	res <- .C("Cdate2jul2",year=as.integer(y),month=as.integer(m),day=as.integer(d),doy=integer(1),PACKAGE="pheno")

	return(list(doy=res$doy, year=res$year))
}
