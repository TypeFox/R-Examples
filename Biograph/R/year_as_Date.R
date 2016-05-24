year_as_Date <-
function(x,format.out)
{   if(missing(format.out)) format.out <- "%Y-%m-%d"  # ISO
	if (!is.numeric(x[1])) stop(paste("year.as.Date: date is not numeric. The (first) date given is: ",x[1],sep=""))
 year <- floor(x)
 k <- as.Date(paste(year,"-01-01",sep=""))
 m <- as.Date(paste(year+1,"-01-01",sep=""))
 z <- as.numeric(k) + (x-year)*as.double(difftime(m,k))+ 10e-10    #10e-12 

 #  as.numeric(as.Date("1970-01-01"))   = 0
 dd <- as.Date(z,origin="1970-01-01")
 #month <- trunc((x-year)*12)+1
 #day <- round((x-year-(month-1)/12)*30.437*12)+1
  date <- format (dd,format.out)
return (date)

}
