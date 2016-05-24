# Number of days between date1 and date2, date as string DD.MM.YYYY.
daysbetween <- function(d1,d2) {
	
	res <- .C("Cdaysbetween",date1=as.character(d1),date2=as.character(d2),ndays=integer(1),PACKAGE="pheno")

	return(res$ndays)
}
