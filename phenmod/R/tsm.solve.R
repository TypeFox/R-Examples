tsm.solve <- function(params, data, silent=FALSE, out2File=FALSE){
	# *****************************************
	return.by.failure <- NA

	# check if data is valid
	data.names <- names(data)
	if (length(which(data.names=="year.lc")) == 0){
		if (!silent){
			cat("Data has to contain column \"year.lc\"!\n")
		}
		return(return.by.failure)
	}
	if (length(which(data.names=="doy.lc")) == 0){
		if (!silent){
			cat("Data has to contain column \"doy.lc\"!\n")
		}
		return(return.by.failure)
	}
	if (length(which(data.names=="latitude")) == 0){
		if (!silent){
			cat("Data has to contain column \"latitude\"!\n")
		}
		return(return.by.failure)
	}
	for (i in 1:365){
		if (length(which(data.names==paste("temperature.",i,sep=""))) == 0){
			if (!silent) {
				cat("Data has to contain column \"temperature.",i,"\"!\n",sep="")
			}
			return(return.by.failure)
		}
	}
	
	# create params as list
	if (!is.list(params)){
		tmp <- params
		params <- list(Tb=tmp[1], Fstar=tmp[2])
	}
	# *****************************************
	
	left.bound <- which(data.names=="temperature.1")
	right.bound <- which(data.names=="temperature.365")

	msg <- ""	

	# PIM on dataframe
	## preparation
	dataframe.length <- dim(data)[1]
	check.factor <- 0.01
	## iteration
	for (dataset.nr in 1:dataframe.length){
		bb.doy.tsm <- 0
		temperatures <- as.numeric(data[dataset.nr, left.bound:right.bound])
		if (is.na(data$doy.lc[dataset.nr])){
			data$doy.bb.tsm[dataset.nr] <- NA
		} else {
			bb.doy.tsm <- .C("tsModel", 
				params=as.numeric(params), 
				lcdoy=as.integer(round(data$doy.lc[dataset.nr])),
				year=as.integer((data$year.lc[dataset.nr]+1)), 
				temperaturevec=temperatures, budburstdoy=as.integer(bb.doy.tsm), 
				PACKAGE="phenmod")$budburstdoy
			if (bb.doy.tsm == 0){
				data$doy.bb.tsm[dataset.nr] <- NA
			} else {
				data$doy.bb.tsm[dataset.nr] <- bb.doy.tsm
			}
		}
		if (!silent){
			if (out2File){
				if (dataset.nr  > (dataframe.length * check.factor)){
					cat(check.factor * 100, "% done!\n",sep="")
					check.factor <- check.factor + 0.01
				}
			} else {
				if(is.na(data$doy.bb.tsm[dataset.nr])){
					cat("\nNo Budburst Doy Found!\n",sep="")
				} else {
					cat(rep("\b", nchar(msg)),sep="")
				}
				msg <- paste(dataset.nr," of ", 
					dataframe.length, " Datapoints done!",sep="")
				cat(msg,sep="")
			}
		}
	}
	if (!silent){
		cat("\n")	
	}
	return(data)
}		
