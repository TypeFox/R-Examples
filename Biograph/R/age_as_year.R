age_as_year <-
function (x,born,format.born)
 { 	if (is.character(x)) x<- as.numeric(x)
	if (!is.numeric(x))
 	   { print ("WARNING in age.as.year: age is not numeric")
 	   	 if (is.character(x)) x <- as.numeric(x) else print ("ERROR in age.as.year: age not numeric and not character.")
 	   }
 	if (missing(format.born)) stop ("age_as_year: format.born is missing")
 	if (missing(born)) stop("ERROR: age.as.year: birth date is missing")
 	if (is.numeric(born)|is.integer(born))
 	  { if (born[1] == 0) {year <- x; return(year)} 
 	  	if (format.born=="CMC") b <- cmc_as_Date (born,1,format.out="%Y-%m-%d")	
 	  	if (format.born=="year"|format.born=="Year"|format.born=="YEAR") d <- born 	
 	  }
 	if (is.character(born)) b <-as.Date(born,format.born)
 	if (class(born)=="Date") b <-as.Date(born,format.born)
    # if (!exists("b")) stop("Error in age_as_year: date of birth cannot be interpreted")
 	if (exists("b"))  d <- Date_as_year (b,format.in="%Y-%m-%d")  # convert date of birth to year
 	year <- d+x
    return (year) 
}
