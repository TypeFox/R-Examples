age_as_Date <-
function (x,born,format.born,format.out)
 { 	if (is.character(x)) x<- as.numeric(x)
	if (!is.numeric(x))
 	   { print ("WARNING in age.as.year: age is not numeric")
 	   	 if (is.character(x)) x <- as.numeric(x) else print ("ERROR in age.as.year: age not numeric and not character.")
 	   }
 	if (missing(format.born)) format.born <- "%Y-%m-%d" # ISO 
 	if (missing(format.out)) format.out <- "%Y-%m-%d" # ISO  
 	if (missing(born)) stop("ERROR: age.as.year: birth date is missing")
 	if (is.character(born)) b <-as.Date(born,format.born)
 	if (class(born)=="Date") b <-as.Date(born,format.born)
 	if (format.born=="CMC"|format.born=="cmc") b<- cmc_as_Date (born,selectday=1,format.out="%Y-%m-%d")
 	d <- Date_as_year (b,format.in="%Y=%m-%d")  # convert date of birth to year
 	year <- d+x
 	date <- year_as_Date (year,format.out=format.out) 
    return (date) 
}
