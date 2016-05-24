cmc_as_Date <-
function (x,selectday,format.out)
 { 	if (is.logical(x)) return(x)  # all values are NA
 	if (!is.numeric(x))
 	   { print ("WARNING in cmc.as.Date: x is not numeric")
 	   	 if (is.character(x)) x <- as.numeric(x) else print ("ERROR in cmc_as_Date: x not numeric and not character.")
 	   }
 	   	
 	if (missing(selectday)) selectday <- 1
 	if (missing(format.out)) format.out  <- "%Y-%m-%d"
 	cmc <- x
 	#cmc[cmc<=0] <- NA
 	# base <- 1900
 	base <- 1000
 	cmc <- 900*12+x
 	year <- base + trunc((cmc-1)/12)
    month <- cmc-(year-base)*12
    day <- rep(selectday,length(cmc))
    year
    month
    z<- ifelse (is.na(cmc),NA,paste(year,"-",month,"-",day,sep=""))
    dd <- as.Date(z,format="%Y-%m-%d")  
    date <- format (dd,format.out)
    return (date)  # format produces a character vector
}
