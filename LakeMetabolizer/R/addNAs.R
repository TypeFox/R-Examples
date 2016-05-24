#'
#'
#'
#'@importFrom utils flush.console
#'@importFrom stats approx
#'
addNAs <- function(x, ...){
	dateL <- grepl(".?date?.", names(x), ignore.case=TRUE) # matches anything with "date" in it, regardless of what else is or is not there
	dL <- grepl("^[dD][oO][yY]$", names(x)) # matches doy, regardless of case
	yL <- grepl("^[yY]ear4?$", names(x))# matches Year, year, year4, Year4

	# =============================================================
	# = The datetime checks result in error if conditions not met =
	# =============================================================
	if(any(dateL)){ # look for the date column
		names(x)[dateL] <- "datetime"
	}else{
		# warning("No 'date' column found")
		stop("No 'datetime' column found")
	}
	if(!"POSIXct"%in%class(x[,"datetime"])){ # make sure the date column is POSIXct (note that if date column is not found, you get this error)
		stop("date column must be POSIXct")
	}
	
	# ===============================================
	# = If doy/ year aren't found, they'll be added =
	# ===============================================
	# Because we are requiring POSIXct datetime, these values can be generated if they're missing
	if(any(dL)){ # look for "day of year column"
		names(x)[dL] <- "doy"
	}else{
		# NOTE: if a "doy" column is not created, byeShort will not work.
		# NOTE: if a "doy" column is not created here, it will have to be created in byeShort anyway
		x[,"doy"] <- date2doy(x[,"datetime"])
		# warning("No 'doy' column found")
		
	}
	if(any(yL)){ # look for "year" column
		names(x)[yL] <- "year"
	}else{
		# NOTE: if a "year" column is not created, byeShort will not work.
		# NOTE: if a "year" column is not created here, it will have to be created in byeShort anyway
		x[,"year"] <- as.integer(format.Date(x[,"datetime"], format="%Y"))
		# warning("No 'year' column found")
	}

	rdig <- 4
	Mode <- function(x){ # note that this function is now in the helper.functions.R file
			ux <- unique(x)
			ux[which.max(tabulate(match(x, ux)))]
	}
	ex <- round(Mode(1/diff(x[,"doy"]))) # note that this is the freq that's calculated in metab.xx()
	mins <- 1/ex*24*60
	is.wholenumber <- function(x, tol = .Machine$double.eps^0.5){abs(x - round(x)) < tol}
	if(!is.wholenumber(mins)){warning("Time between samples not whole number")}
	x1 <- byeShort(X=x, Expected=ex, ToCount="doy", TruncToCount=TRUE, ...)
	if(nrow(x1)==0){
		return(x1)
	}
	Range <- range(x1[,"datetime"]) #intentionally not truncating (having already used byeShort, I don't have to worry about starting and stopping at a specific time each day)
	Ideal <- data.frame("datetime"=seq(Range[1], Range[2], by=paste(mins, "mins")))
	
	print(paste("NA's added to fill in time series:",dim(Ideal)[1]-dim(x1)[1], sep=" "))
	flush.console()
	x2 <- merge(x1, Ideal, all.x=TRUE, all.y=TRUE)
	if(any(yL)){x2[,"year"] <- approx(x2[,"datetime"], x2[,"year"], xout=Ideal[,1])$y}
	if(any(dL)){x2[,"doy"] <- approx(x2[,"datetime"], x2[,"doy"], xout=Ideal[,1])$y}
	
	return(x2)
}
