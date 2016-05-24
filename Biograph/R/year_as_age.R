year_as_age <-
function (x,born,format.born)
 { 	# Decimal year
	if (is.character(x)) x<- as.numeric(x)
	if (!is.numeric(x))
 	   { print ("WARNING in year.as.age: year is not numeric")
 	   	 if (is.character(x)) x <- as.numeric(x) else print ("ERROR in year.as.age: year not numeric and not character.")
 	   }
 	if (missing(format.born)) format.born <- "%Y-%m-%d" # ISO 
 	if (missing(born)) stop("ERROR: year.as.age: birth date is missing")
 	if (class(born)=="Date" | is.character(born)) b <-as.Date(born,format.born)
 	if (format.born=="year" | is.numeric(born))
 	   { age = (x-born)
 	   	 return (age)
 	   }
 	if (is.numeric(born))
 	  { age <- x - born }   else
 	  { d <- Date_as_year (born,format.in=format.born)  # convert date of birth to year
 	    age <- x - d } 
 	k <- age<0
 	if (TRUE%in%k) warning ("At least one date (year) is before date of birth.",call.=TRUE)
    return (age) 
}
