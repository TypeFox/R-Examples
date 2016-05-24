ddate <-
function(year, month, day)
{
    # convert year, month, day to decimal date
    if(length(year) > 1) {
	month <- year[2]
	day <- year[3]
	year <- year[1]
    }
    year + julian(as.Date(paste(year, month, day, sep="-")), origin=as.Date(paste(year-1, "12", "31", sep="-")))[[1]]/365
    
}

