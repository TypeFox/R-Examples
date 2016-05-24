#################################################################################
##
##   R package parma
##   Alexios Ghalanos Copyright (C) 2012-2013 (<=Aug)
##   Alexios Ghalanos and Bernhard Pfaff Copyright (C) 2013- (>Aug)
##   This file is part of the R package parma.
##
##   The R package parma is free software: you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
##   the Free Software Foundation, either version 3 of the License, or
##   (at your option) any later version.
##
##   The R package parma is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##
#################################################################################



# generate non weekend days given a start and finish
SeqWeekdays = function(start, end, format = "%d-%m-%Y"){
	f = format(strptime(as.character(start), format), format = "%Y-%m-%d")
	e = format(strptime(as.character(end), format), format = "%Y-%m-%d")	
	z = seq.Date(from = as.Date(f), to = as.Date(e), by = "days", origin. = "1900-01-01")
	z = z[-c(which(weekdays(z) == "Saturday"))]
	z = z[-c(which(weekdays(z) == "Sunday"))]
	z = as.character(z)
	attr(z, "format") = "%Y-%m-%d"
	return( z )
}

# choose a weekday from a set of dates
Weekday = function(dates, format = "%Y-%m-%d", weekday = "Friday"){
	z = as.Date(dates, format = format)
	z = z[which(weekdays(z) == weekday)]
	z = as.character(z)
	attr(z, "format") = "%Y-%m-%d"
	return( z )
}

LastMonthDay = function(dates, format = "%Y-%m-%d"){
	z = as.Date(dates, format = format)
	idx = which(diff(as.integer( format(z, "%m") ) )!=0)
	return( idx )
}

Lag = function(x, n.lag = 1, pad = 0){
	if(is.numeric(x)){
		return(.numeric.lag(x, n.lag, pad))
	} else if(is.matrix(x)){
		return(.matrix.lag(x, n.lag, pad))
	} else if(is(x, "timeSeries")){
		return(.timeSeries.lag(x, n.lag, pad))
	} else{
		stop("\nunrecognized/unsupported format for x.")
	}
}

.numeric.lag = function(x, n.lag = 1, pad = 0)
{
	if(length(n.lag) == 1){
		ans = c(rep(pad, n.lag), head(x, length(x) - n.lag))
	} else{
		ans = NULL
		ans = apply(as.data.frame(n.lag), 1, FUN = function(i) rbind(ans, c(rep(pad, i), head(x, length(x) - i))))
	}
	return(ans)
}

.matrix.lag = function(x, n.lag = 1, pad = 0)
{
	if(length(n.lag) == 1){
		ans = rbind(matrix(pad, ncol = dim(x)[2], nrow = i), head(x, dim(x)[1] - n.lag))
	} else{
		ans = NULL
		for(i in n.lag){
			ans = cbind(ans, rbind(matrix(pad, ncol = dim(x)[2], nrow = i), head(x, dim(x)[1] - i)))
		}
	}
	return(ans)
}

.timeSeries.lag = function(x, n.lag = 1, pad = 0)
{
	cnames = colnames(x)
	if(length(n.lag) == 1){
		cxnames = paste(cxnames, ".lag", n.lag, sep = "")
		ans = rbind(matrix(pad, ncol = dim(x)[2], nrow = i), head(x, dim(x)[1] - n.lag))
	} else{
		cxnames = NULL
		for(i in n.lag) cxnames = c(cxnames, paste(cnames, ".lag", i, sep = ""))
		ans = NULL
		for(i in n.lag){
			ans = cbind(ans, rbind(matrix(pad, ncol = dim(x)[2], nrow = i), head(x, dim(x)[1] - i)))
		}
	}
	return(ans)
}