# Converts a string date "DD.MM.YYYY" into a Julian date (doy,year)
date2jul1 <- function(d) {
	
	res <- .C("Cdate2jul1",date=as.character(d),doy=integer(1),year=integer(1),PACKAGE="pheno")

	return(list(doy=res$doy, year=res$year))
}
