#################################################################################
##
##   R package rugarch by Alexios Ghalanos Copyright (C) 2008-2015.
##   This file is part of the R package rugarch.
##
##   The R package rugarch is free software: you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
##   the Free Software Foundation, either version 3 of the License, or
##   (at your option) any later version.
##
##   The R package rugarch is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##
#################################################################################

################################################################################
# EXPORTED
move = function(index, by=1){
	if(!is(index, "POSIXct") && !is(index, "Date") && !is(index, "numeric"))
		stop("\nunrecongnized time index")
	n = length(index)
	if(by>=n) stop("\nby must be less than length of index!")
	if(by==0) return(index)
	p = median(diff(index))
	newindex = c(index[-c(1:by)], generatefwd(index[n], length.out=1, by=p))
	return(newindex)
}

generatefwd = function(T0, length.out = 1, by = "days"){
	if(!is(T0, "POSIXct") && !is(T0, "Date") && !is(T0, "numeric"))
		stop("\nunrecongnized time index")
	Z = seq(T0, by = by, length.out=length.out*4)[-1]
	if(!is.numeric(T0)){
		W = weekdays(Z)
		idx = c(which(W=="Saturday"), which(W=="Sunday"))
		if(length(idx)>0) Z = Z[-idx]
	}
	Z = Z[1:length.out]
	return(Z)
}

# T0 = as.POSIXct("2001-01-01 16:00:00")
# interval = format(seq(as.POSIXct("2001-01-01 09:30:00"), as.POSIXct("2001-01-01 16:00:00"), by="min"), "%H:%M:%S")
# by = "mins"
# length.out=1000
ftseq = function(T0, length.out, by, interval, exclude.weekends = TRUE)
{
	start = T0
	# Just one check:
	if(!is(start, "POSIXct")) stop("\nstart must be a POSIXct object")
	U = format(start, "%H:%M:%S")
	if(is.na(match(U, interval))) stop("\nstart must match one of the supplied interval values.")
	zn = length.out-1
	k = 1
	while(zn<length.out){
		# function does not know what increment "by" is nor the interval.
		# start out with length.out*10 as an estimate and continue looping
		# until we have enough points to satisfy the requirements, including
		# exclusion of weekends
		z = seq(from = start, by = by, length.out = length.out*10*k)[-1]
		y = format(z, "%H:%M:%S")
		z = z[which(!is.na(match(y, interval)))]
		if(exclude.weekends){	
			z = z[-which(weekdays(z)=="Saturday")]
			z = z[-which(weekdays(z)=="Sunday")]
		}
		zn = length(z)
		k = k*2
	}
	return(z[1:length.out])
}

#R = ftseq(as.POSIXct("2001-01-01 16:00:00"), length.out=2000, by = 18000, interval = interval)
################################################################################
.extractdata = function(data, warn = FALSE)
{
	xdata = try(as.xts(data), silent = TRUE)
	if(inherits(xdata, "try-error")){
		if(warn) warning("\nrugarch-->warning: data indexing not recognized by xts...coercing to Date with origin 1970-01-01.")
		if(is.data.frame(data) | is.matrix(data)) data = as.numeric(data[,1]) else data = as.numeric(data)
		data = unname(data)
		xdata = xts(data, as.POSIXct(as.Date(seq_along(data), origin="1970-01-01")))
	}
	obj = list()
	obj$data = as.numeric(coredata(xdata))
	obj$index = index(xdata)
	obj$period = median(diff(index(xdata)))
	return(obj)
}

.numeric2xts = function(data){
	data = as.numeric(data)
	return(xts(data, as.Date(1:NROW(data), origin="1950-01-01")))
}

.genxts = function(index0, length.out = 10, period = "days"){
	Z = seq(index0, by = period, length.out=length.out)
	return(Z)
}
################################################################################
# DEPRECATED
# make generatefwd search for available dates from original set
.generatefwd = function(Dates, N, dformat, periodicity = "days")
{
	if(is.numeric(Dates[1])){
		n = length(Dates)
		fwd = (Dates[n]+1):(n+N)
	} else{
		n = length(Dates)
		# reformat data
		sdx = format(Dates[n], "%m/%d/%y")
		fwd = chron::seq.dates(from=sdx, by = periodicity, length.=N*8)
		# generate enough to dates to not cause problem with weekend exclusion later
		z1 = which(chron::is.weekend(fwd))
		fwd = as.character(fwd)
		fwd = fwd[-z1]
		fwd = fwd[2:(N+1)]
		fwd = as.character(format(strptime(as.character(fwd), "%m/%d/%y"), dformat))
		fwd = as.Date(fwd, format = dformat)
	}
	return(fwd)
}

#ForwardDates = function(Dates, n.ahead, date.format, periodicity="days")
#{#
#	UseMethod("ForwardDates")
#}

.ForwardDates = function(Dates, n.ahead, date.format, periodicity="days")
{
	if(is.numeric(Dates[1])){
		n = length(Dates)
		fwd = (Dates[n]+1):(n+n.ahead)
	} else{
		n = length(Dates)
		# reformat data
		sdx = format(as.Date(Dates[n], format = date.format), format = "%m/%d/%y")
		fwd = chron::seq.dates(from = sdx, by = periodicity, length. = 8 * n.ahead)
		# generate enough to dates to not cause problem with weekend exclusion later
		z1 = which(chron::is.weekend(fwd))
		fwd = format(as.Date(as.character(fwd), format = "%m/%d/%y" ), format = date.format)
		if(length(z1)>0) fwd = fwd[-z1]
		fwd = fwd[2:(n.ahead+1)]
		fwd = as.Date(fwd, format = date.format)
	}
	return(fwd)
}
	
#setMethod("ForwardDates", signature(Dates = "character"), definition=.ForwardDates)

#WeekDayDummy = function(Dates, date.format, weekday = "Monday")
#{
#	UseMethod("WeekDayDummy")
#}

.WeekDayDummy = function(Dates, date.format, weekday = "Monday")
{
	valid.weekdays = c("Monday","Tuesday","Wednesday","Thursday","Friday")
	wd = match(weekday, valid.weekdays)
	if(is.na(wd)) stop("not a valid weekday")
	xdates = .makedate(Dates)
	if(xdates$status == 0) stop("\nForwardDates-->error: date format not recognized\n", call. = FALSE)
	ddates = xdates$dates
	dformat = xdates$dformat
	ans = rep(0, length(Dates))
	zz = which(weekdays(ddates) == weekday)
	ans[zz]=1
	ans
}

#setMethod(f = "WeekDayDummy", signature(Dates = "character"), definition = .WeekDayDummy)

# create the forecast date set which depends on out.sample data
.forcdates = function( origdates, n.ahead, N, i, ns , dformat)
{
	# need to return both a timevector and character vector
	if( ns == 0 ){
		fwwdx = .generatefwd(origdates, N = n.ahead, dformat = dformat, 
				periodicity = "days")
	} else{
		if( (n.ahead + i - 1) <= ns ){
			fwwdx = origdates[(N + i):(N + i + n.ahead -1)]
		} else{
			if( i <= ns ){
				nt = max(0, ns - i + 1)
				fwwdx1 = origdates[(N+i):(N+ns)]
				fwwdx2 = .generatefwd(origdates[1:(N+ns)], N = n.ahead-nt, 
								dformat = dformat, periodicity = "days")
				fwwdx =c(fwwdx1, fwwdx2)		
			} else{
				fwwdx = .generatefwd(origdates[1:(N+ns)], N = n.ahead, 
								dformat = dformat, periodicity = "days")
			}
			
		}
	}
	fwdd = fwwdx
	return(fwdd)
}


.makedate = function(x)
{
	# find the divisor: 4 cases "-", "/", ".", and no divisor
	allc = strsplit(x[1], "")
	
	if(any(allc[[1]] == "-")){
		dt = "-"
		dte = t(apply(as.data.frame(x), 1, FUN=function(z) as.numeric(strsplit(z, dt)[[1]]) ))
		ld = max(apply(dte, 1, FUN = function(x) max(nchar(as.character(x)))))+3		
	} else if(any(allc[[1]] == "/")){
		dt = "/"
		dte = t(apply(as.data.frame(x), 1, FUN=function(z) as.numeric(strsplit(z, dt)[[1]]) ))
		ld = max(apply(dte, 1, FUN = function(x) max(nchar(as.character(x)))))+3
	} else if(any(allc[[1]] == ".")){
		dt = "."
		dte = t(apply(as.data.frame(x), 1, FUN=function(z) as.numeric(strsplit(z, dt)[[1]]) ))
	} else{
		# this is a little more complicated
		ld = length(allc[[1]])
		if(ld==6){
			dte = t(apply(as.data.frame(x), 1, FUN=function(z) 
								as.numeric(c(substr(z, 1,2), substr(z, 3,4), substr(z, 5,6)))))
		} else if(ld==8){
			# 2 cases either the 4 digit year is at the beginning or else at the end
			dte.1 = as.vector(t(apply(as.data.frame(x), 1, FUN=function(z) 
								as.numeric(c(substr(z, 1,2))))))
			dte.2 = as.vector(t(apply(as.data.frame(x), 1, FUN=function(z) 
										as.numeric(c(substr(z, 5,6))))))
			if(all(dte.1>18)){
				dte = t(apply(as.data.frame(x), 1, FUN=function(z) 
									as.numeric(c(substr(z, 1,4), substr(z, 5,6), substr(z, 7,8)))))
			} else if(all(dte.2>18)){
				dte = t(apply(as.data.frame(x), 1, FUN=function(z) 
									as.numeric(c(substr(z, 1,2), substr(z, 3,4), substr(z, 5,8)))))
			} else{
				return(list(status=0))
			}
		} else{
			return(list(status=0))	
		}
	}
	m = 0
	for(i in 1:3){
		if(all(dte[,i]<=12)) m = i
	}
	if(m==0) return(list(status=0))
	sq = 1:3
	sq = sq[-m]
	y = 0 
	for(i in sq){
		if(any(dte[,i]>31)) y = i
	}
	if(y==0) return(list(status=0))
	d = sq[-which(sq==y)]
	dmatrix = cbind(dte[,d], dte[,m], dte[,y])
	colnames(dmatrix) = c("d","m","y")
	if(ld==6){
		ddates = as.Date(paste(dmatrix[,3], dmatrix[,2], dmatrix[,1], sep = "-"), format="%y-%m-%d")
		dformat = "%y-%m-%d"
	} else{
		ddates = as.Date(paste(dmatrix[,3], dmatrix[,2], dmatrix[,1], sep = "-"), format="%Y-%m-%d")
		dformat = "%Y-%m-%d"
	}
	

	return(list(datesmat = dmatrix, dates = ddates, dformat = dformat, status=1))
}