Date_as_year <-
function (x,format.in)
 {  if (missing(format.in)) format.in <- "%Y-%m-%d"
 	if (substr(format.in,1,1)!="%") stop ("Date_as_year: format.in is not date format")
 	if (is.character(x)) x<- as.Date(x,format.in)
 	# Check whether date is feasible
 	d2 <- as.POSIXlt(x)
 	year <- d2$year+1900
 	# MicMac year <- as.numeric(unlist(strsplit(as.character(d2),"-"))[1])
 	# see MicPostFun function diffBetweenDates
    if (length(na.omit(year))==0) {year.frac <- rep(NA,length(year)) ; return(year.frac)}
    if (min(year,na.rm=TRUE) < 32) stop(paste("Date_as_year: Please check date format. Year = ",year,sep=""))
 	# Julian dates
 	k <- ifelse (is.na(year),NA,paste(year,"-01-01",sep=""))
 	k <- as.Date(k)
 	m <- ifelse (is.na(year),NA,paste(year+1,"-01-01",sep=""))
 	m <- as.Date(m)
 	# Fraction of year
 	frac <- (as.double(difftime(x,k,units="days")))/(as.double(difftime(m,k,units="days"))) # Diff in days
 	year.frac <- year + frac
 	

    #year.frac <- year+(month-1)/12+(day-1)/(30.437*12)  # check see hulp.r
    return (year.frac)  # date is object of class 'Date'
}
