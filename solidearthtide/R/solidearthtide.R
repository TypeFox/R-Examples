#	*	*	*	Function "CalculateTide" Package "solidearthtide"		*	*	*
CalculateTide <- function(fname, yearIn, monthIn, dayIn, latIn, lonIn, boolRet=FALSE)
{
if (!is.character(fname)) stop('Argument <<fname>> should be a string')
if (!is.numeric(yearIn)) stop('Argument <<yearIn>> should be a number')
if (!is.numeric(monthIn)) stop('Argument <<monthIn>> should be a number')
if (!is.numeric(latIn)) stop('Argument <<latIn>> should be a number')

# iyr	year    [1980-2016]
# imo	month number [1-12]
# idy	day          [1-31]
# glad	Lat. (pos N.) [- 90, +90]
# glod	Lon. (pos E.) [-360,+360]
retArray <- matrix(0,1441, 4)
retErr <- 0
ret <- .Fortran("calctide",iyr = as.integer(yearIn), imo = as.integer(monthIn), 
idy = as.integer(dayIn), glad = as.double(latIn), glod = as.double(lonIn), iretarray=as.integer(boolRet), retArray=retArray,
retErr=as.integer(retErr) ,DUP = TRUE, PACKAGE="solidearthtide")
if (retErr !=0) stop(paste('Error ',retErr))
if (fname != '') write.table(ret$retArray, fname, quote=FALSE, sep=" ", row.names=FALSE) # save the output array to a file
if (boolRet) return(ret$retArray) # return the output array
}

Datenum <- function(d) 
{# Convert date to serial date number with the same values as Matlab's datenum
return (as.integer(719529+(d-as.Date("1970-01-01"))))
}

Datevec <- function(s)
{# Convert serial date number to date with the same values as Matlab's datevec
(s-719529)+as.Date("1970-01-01")
}

DatevecV <- function(s)
{# Convert serial date number to list(year,month,day) with the same values as Matlab's datevec
dTMP<-(s-719529)+as.Date("1970-01-01")
list(year=as.POSIXlt(dTMP)$year + 1900,month=as.POSIXlt(dTMP)$mon + 1,day=as.POSIXlt(dTMP)$mday)
}

GetLeapDates <- function()
{# return the leap dates from the UTC to GPS Time Conversion Table
data("UTCtoGPSTimeConversion", envir = environment())
UTCtoGPSTimeConversion<-get("UTCtoGPSTimeConversion", envir  = environment())
leapDates <- UTCtoGPSTimeConversion[which(!is.na(UTCtoGPSTimeConversion[["GPSminusUTC"]])), "yearLimitFrom"]
leapDates[1] <- "1980-01-06"
return(leapDates)
}

UTC2GPS <- function(utctime, asMATLAB=FALSE)
{
# Convert UTC(GMT) time tags to GPS time accounting for leap seconds
# Based on the MATLAB utc2gps Copyright 2008 Ian M. Howat
# this version can return MATLAB serial date number or R date object
stepdates <- GetLeapDates()# get the leap dates
stepdates <- Datenum(as.Date(stepdates, "%Y-%m-%d"))
steptime <- 0:(length(stepdates)-1)/86400
utctime <- kronecker(matrix(1,1,length(stepdates)),utctime)
stepdates<-kronecker(matrix(1,dim(utctime)[1],1),matrix(stepdates,1))
tmp <- utctime - stepdates
tmp <- ifelse(tmp>0.2,1,0)
whStepTime<- apply(tmp,1,sum)
gpstime<-utctime[,1]  + steptime[whStepTime]
if (asMATLAB) return(gpstime) else return(as.Date(gpstime-719529))
}

Earthtide <- function(latIn, lonIn, timeSeq, boolHorizontalMotion=FALSE)
{# predict displacements from solid earth tides
# based on 
uniqueDays<-unique(floor(timeSeq))
upOut<-matrix(0,1,length(floor(timeSeq)))
if (boolHorizontalMotion){
northOut<-upOut
eastOut<-upOut
}
for (k in 1:length(uniqueDays))
{
indUniq <- ifelse(floor(timeSeq)==uniqueDays[k],1,0)
vDate<-DatevecV(uniqueDays[k])
tideData<-CalculateTide('', vDate[["year"]], vDate[["month"]], vDate[["day"]], latIn, lonIn, TRUE)
ups = tideData[,4] #OK
ts <-Datenum(as.Date(paste(vDate[["year"]], vDate[["month"]], vDate[["day"]],sep=' '), "%Y %m %d"))
ts <-ts+tideData[,1]/86400
#print(format(round(ts, 7), nsmall = 7)) #OK
inU<-c(indUniq,rep(0,1441-337))
upOut[1,which(indUniq==1)]<-spline(ts,ups,sum(indUniq))$y
if (boolHorizontalMotion) {
ups <- tideData[,2]
northOut[1,which(indUniq==1)]<-spline(ts,ups,sum(indUniq))$y
ups <- tideData[,3]
eastOut[1,which(indUniq==1)]<-spline(ts,ups,sum(indUniq))$y
}
}
if (boolHorizontalMotion) return (list(upOut=upOut, northOut=northOut, eastOut=eastOut)) else return(upOut)
}


