# Author: cns
# advanced models

dlmodeler.countweekdays <- function(days)
{
	days <- weekdays(days)
	count.fn <- function(wd) sum(days==wd)
	wds <- list("Monday","Tuesday","Wednesday","Thursday","Friday","Saturday","Sunday")
	lapply(wds,count.fn)
}

dlmodeler.yeardays <- function(year)
{
	days <- seq(as.Date(paste(year,"01-01",sep="-")),as.Date(paste(year+1,"01-01",sep="-"))-1,by="day")
	return(dlmodeler.countweekdays(days))
}
