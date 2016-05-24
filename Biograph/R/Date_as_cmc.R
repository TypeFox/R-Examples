Date_as_cmc <-
function (x,format.in)
  { if (missing(format.in)) format.in <- "%Y-%m-%d"
  	if (is.character(x)) d <- as.Date (x,format.in)
  	if (class(x) =="Date") d <- x
  	lt <- as.POSIXlt(as.Date(d, origin="1900-01-01"))
	cmc <- lt$year*12 + lt$mon + 1
	return(list (cmc=cmc,
	             selectday=lt$mday))
  }
