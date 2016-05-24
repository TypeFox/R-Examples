Date_as_age <-
function (x,format.in,born)
 {  if (missing(format.in)) stop("Date_as_age: format.in is missing")
 	if (missing (born)) stop("Date_as_age: date of birth is missing")
  	if (substr(format.in,1,1)!="%") stop ("Date_as_age: format.in is not date format")
 	if (is.character(x)) x<- as.Date(x,format.in)
 	d <- as.POSIXlt(x,format=format.in,tz="UTC")
 	d.born <- as.POSIXlt (born,format=format.in,tz="UTC")
 	age.days=(d-d.born)  # age in days 	
 	
 #	require (lubridate)
 	z <- as.interval (age.days,d.born)
    p=as.period(z,unit="year")
    p@year
    p@month
    p@day
    age = as.duration (z)
    age.sec <- age@.Data  # age in seconds
    days = age.sec/86400  # = age.days
    year = Date_as_year(d,"%Y-%m-%d")
    birth = Date_as_year(d.born,"%Y-%m-%d")     
    year.average = days/(year-birth)  # Average length of year in days during this period
     
    return (list (age.sec=age.sec,
                  age.days=days,
                  age = p,
                  age.year=year-birth )) 
}
