pim.solve <- function(params, data, model.no=1, silent=FALSE, out2File=FALSE){
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

	# check if model.no is valid
	if ((model.no < 1)||(model.no > 12)){
		if (!silent){
			cat("Invalid model number! (1 <= model.no <= 12)\n",sep="")
		}
		return(return.by.failure)
	}
	
	# create params as list
	if (!is.list(params)){
		tmp <- params
		params <- list(a1=tmp[1], a2=tmp[2], a3=tmp[3], a4=tmp[4], 
			T.min.i=tmp[5], T.opt.i=tmp[6], T.max.i=tmp[7],
			T.min.p=tmp[8], T.opt.p=tmp[9], T.max.p=tmp[10])
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
		bb.doy.pim <- 0
		temperatures <- as.numeric(data[dataset.nr, left.bound:right.bound])
		if (is.na(data$doy.lc[dataset.nr])){
			data$doy.bb.pim[dataset.nr] <- NA
		} else {
			bb.doy.pim <- .C("pimModel", 
				params=as.numeric(params), 
				lcdoy=as.integer(round(data$doy.lc[dataset.nr])),
				year=as.integer((data$year.lc[dataset.nr]+1)), 
				latitude=as.numeric(data$latitude[dataset.nr]),
				temperaturevec=temperatures,modelno=as.integer(model.no), 
				budburstdoy=as.integer(bb.doy.pim), 
				PACKAGE="phenmod")$budburstdoy
			if (bb.doy.pim == 0){
				data$doy.bb.pim[dataset.nr] <- NA
			} else {
				data$doy.bb.pim[dataset.nr] <- bb.doy.pim
			}
		}
		if (!silent){
			if (out2File){
				if (dataset.nr  > (dataframe.length * check.factor)){
					cat(check.factor * 100, "% done!\n",sep="")
					check.factor <- check.factor + 0.01
				}
			} else {
				if(is.na(data$doy.bb.pim[dataset.nr])){
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