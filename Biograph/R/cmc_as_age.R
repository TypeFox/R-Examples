cmc_as_age <-
function (x,born,format.born)
 { 	# Decimal year
	if (is.character(x)) x<- as.numeric(x)
	if (!is.numeric(x))
 	   { print ("WARNING in cmc.as.age: cmc is not numeric")
 	   	 if (is.character(x)) x <- as.numeric(x) else print ("ERROR in cmc.as.age: cmc not numeric and not character.")
 	   }
 	if (missing(format.born)) stop ("cmc_as_age: format.born is missing")
 	if (missing(born)) stop("ERROR: cmc_as_age: birth date is missing")
 	if (class(born)=="Date" | is.character(born)) b <-as.Date(born,format.born)
 	if (format.born=="CMC"|format.born=="cmc")
 	   { age = (x-born)/12
 	   	 if (length(which (age < 0 )) > 0) # At least one negative age 
 	   	  { age[age < 0] <- NA 
 	   	  	message ("cmc_as_age: negative ages are replaced by NA")  }
 	   	  	
 	   	 c <- cmc_as_year (x)
 	   	 return (list(year = c,
 	   	              age = age))
 	   }
 	d <- Date_as_year (born,format.in="%Y-%m-%d")  # convert date of birth to year
 	print (d)
 	print (born)
 	c <- cmc_as_year (x)
 	z <- c - d
 	k <- z<0
 	if (TRUE%in%k) warning ("At least one date (cmc) is before date of birth.")
 	age <- c - d
    return (list (year = c,
                  age = age)) 
}
