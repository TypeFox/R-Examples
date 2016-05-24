#'@name metab
#'@aliases 
#'metab
#'@title Calculate metabolism
#'@description 
#'Returns daily time series of gross primary production (GPP), respiration (R), and net ecosystem production (NEP). Depending on the method used, other information may be returned as well. Calculations are made using one of 5 statistical methods.
#'
#'
#'@usage
#'metab(data, method, wtr.name="wtr", irr.name="irr", do.obs.name="do.obs", ...)
#'
#'@param data 
#' a data.frame whose columns are 
#'
#'  \code{"datetime"} = class POSIXct vector
#'
#'  \code{"do.obs"} = numeric vector of oxygen concentration in mg/L
#'
#'  \code{"do.sat"} = numeric vector of saturated oxygen concentration in mg/L
#'
#'  \code{"k.gas"} = numeric vector of gas exchange coefficient values in m/day, should be 0 when depth of do.obs is deeper than z.mix
#'
#'  \code{"z.mix"} = numeric vector of mixing depth values in meters
#'
#'  \code{"irr"} = numeric vector of PAR values, arbitrary units
#'
#'  \code{"wtr"} = numeric vector of water temperature values, arbitrary units
#'
#' Columns that are not used by a particular statistical method do not need to be supplied.
#'
#'@param method
#' a character string specifying one of the 5 statistical methods 
#'(bayesian, bookkeep, kalman, ols, mle)
#'@param wtr.name the name of the column containing temperature at the depth of do.obs (predictor variable for R)
#'@param irr.name the name of the column containing irradiance (predictor variable for GPP)
#'@param do.obs.name the name of the column in data containing the DO observations (in mg/L) to be used as the response variable
#'@param ... arguments to be passed on to the metabolism model specified by \code{method}
#'
#' @return 
#' A data.frame containing columns for year, doy (day of year, julian day plus fraction of day), GPP, R, and NEP
#' \item{year}{integer year}
#' \item{doy}{numeric, day of year + fraction of day, where the day is the julian day, and a fraction of 0.5 corresponds to noon}
#' \item{GPP}{numeric, gross primary production, in units of mg O2 per liter per day. By convention, this value is positive.}
#' \item{R}{numeric, respiration, in units of mg O2 per liter per day. By convention, this value is negative}
#' \item{NEP}{numeric, net ecosystem production, in units of mg O2 per liter per day. For most methods this equal GPP+R, but this is not necessarily the case for \code{"method"="bookkeep"}}
#' Note that different models will have different \link{attributes} attached to them. See examples. 

#'@keywords metabolism
#'
#'@author
#'Ryan D. Batt
#'
#'@seealso 
#' Metabolism models: \link{metab.bookkeep}, \link{metab.ols}, \link{metab.mle}, \link{metab.kalman}, \link{metab.bayesian}
#' 
#' For smoothing noisy temperature: \link{temp.kalman}
#'
#' To calculate do.sat: \link{o2.at.sat}
#'
#' To calculate k.gas: \link{k600.2.kGAS}
#'
#' To calculate k600 values for k.gas: \link{k.cole}, \link{k.crusius}, \link{k.macIntyre}, \link{k.read}
#'
#'
#'@examples
#'
#'	# fake data
#'	datetime <- seq(as.POSIXct("2014-06-16 00:00:00", tz="GMT"), 
#'      as.POSIXct("2014-06-17 23:55:00", tz="GMT"), length.out=288*2)
#'	do.obs <- 2*sin(2*pi*(1/288)*(1:(288*2))+1.1*pi) + 8 + rnorm(288*2, 0, 0.5)
#'	wtr <- 3*sin(2*pi*(1/288)*(1:(288*2))+pi) + 17 + rnorm(288*2, 0, 0.15)
#'	do.sat <- LakeMetabolizer:::o2.at.sat.base(wtr, 960)
#'	irr <- (1500*sin(2*pi*(1/288)*(1:(288*2))+1.5*pi) +650 + rnorm(288*2, 0, 0.25)) * 
#'      ifelse(is.day(datetime, 42.3), 1, 0)
#'	k.gas <- 0.4
#'	z.mix <- 1
#'
#' # plot time series
#'	plot(wtr, type="l", xaxt="n", yaxt="n", xlab="", ylab="")
#'	par(new=TRUE); plot(do.obs, type="l", col="blue", xaxt="n", yaxt="n", xlab="", ylab="")
#'	par(new=TRUE); plot(irr, type="l", col="orange", xaxt="n", yaxt="n", xlab="", ylab="")
#'	abline(v=144, lty="dotted")
#'	abline(v=288)
#'	legend("topleft", legend=c("wtr", "do.obs", "irr"), lty=1, 
#'      col=c("black", "blue", "orange"), inset=c(0.08, 0.01))
#'
#' # put data in a data.frame
#'	data <- data.frame(datetime=datetime, do.obs=do.obs, do.sat=do.sat, k.gas=k.gas, 
#'      z.mix=z.mix, irr=irr, wtr=wtr)
#'
#' # run each metabolism model
#'	m.bk <- metab(data, "bookkeep", lake.lat=42.6)
#'  m.bk <- metab(data, lake.lat=42.6) # no method defaults to "bookeep" 
#'	m.ols <- metab(data, "ols", lake.lat=42.6)
#'	m.mle <- metab(data, "mle", lake.lat=42.6)
#'	m.kal <- metab(data, "kalman", lake.lat=42.6)
#'	\dontrun{m.bay <- metab(data, "bayesian", lake.lat=42.6)}
#' 
#' # example attributes
#' names(attributes(m.ols))
#' attr(m.ols, "mod")
#'  
#' # To get full JAGS model
#' # including posterior draws:
#' \dontrun{names(attributes(m.bay))}
#' \dontrun{attr(m.bay, "model")}
#'
#'@export

metab <- function(data, method = NULL, wtr.name="wtr", irr.name="irr", do.obs.name="do.obs", ...){
	
	m.args <- list(...)
	
	#Rename the WTR column to be used (must be wtr to easily be passed to meta.* functions)
	if(wtr.name != "wtr"){
		if(!"wtr"%in%names(data)){
			names(data)[names(data)==wtr.name] <- "wtr"
		}else{
			data[,"wtr"] <- data[,wtr.name]
		}
	}
	
	if(irr.name != "irr"){
		if(!"irr"%in%names(data)){
			names(data)[names(data)==irr.name] <- "irr"
		}else{
			data[,"irr"] <- data[,irr.name]
		}
	}
	
	if(do.obs.name != "do.obs"){
		if(!"do.obs"%in%names(data)){
			names(data)[names(data)==do.obs.name] <- "do.obs"
		}else{
			data[,"do.obs"] <- data[,do.obs.name]
		}
	}
	
	
	# ===================
	# = Identify method =
	# ===================
	possibleMethods <- c("bookkeep", "bayesian",  "kalman", "ols", "mle")
	mtd <- match.arg(method, possibleMethods)
	
	mtdCall <- paste("metab",mtd, sep=".")
	
	# ==============
	# = Groom data =
	# ==============
	# Removes days with many NA's:
	data1 <- addNAs(data[complete.cases(data),], percentReqd=1) # note that addNAs ALSO checks for POSIXct datetime, and adds year/doy
	data2 <- data1[complete.cases(data1),]
	
	# ==================================
	# = Prepare to apply metab to data =
	# ==================================
	ids <- id(list(data2[,"year"],trunc(data2[,"doy"]))) # ID chunks to be analyzed
	ids <- as.integer(ids - (min(ids)-1))
	nid <- length(unique(ids))
	results <- vector("list", nid)
	
	# ==================================
	# = Apply metab to subsets of data =
	# ==================================
	for(i in unique(ids)){
		
		poss.args <- c("do.obs","do.sat","k.gas","z.mix", "irr", "wtr", "datetime") # data2 columns that could correspond to arguments
		used.args <- poss.args[poss.args%in%names(data2)] # assuming arguments are used if they are in data2
		largs0 <- as.list(data2[i==ids, used.args]) # subsetting columns of data2 to used.args (just a safety check, probably not needed)
		largs <- c(largs0, m.args[!names(m.args)%in%names(largs0)]) # adding on any other arguments supplied via ...
		# note that in largs, argument supplied through data/data2/poss.args take precedent over arguments from ...
		
		# print(paste("Analyzing day #", i)); flush.console(); # Is this annoying? I'm commenting-out
		results[[i]] <- do.call(mtdCall, largs) # this is where all of the work happens
	}
	answer0 <- conquerList(results, naming=data.frame("year"=data2[!duplicated(ids),"year"], "doy"=trunc(data2[!duplicated(ids),"doy"])))
	
	
	a0.names <- names(results[[1]])
	
	# =======================================================
	# = Add non-metab list elements as attributes to output =
	# =======================================================
	# only need to add attributes if it's a list (more than 1 element, not a data frame)
	if(length(a0.names)>1 & is.list(answer0) & !is.data.frame(answer0)){
		
		names(answer0) <- a0.names
		answer <- answer0$metab
		for(i in 1:length(a0.names)){
			if(a0.names[i]=="metab"){next}
			if(a0.names[i]=="smoothDO"){ # do a little extra attribute work if smoothDO
				t.sDO <- answer0[[a0.names[i]]] # grab the element of the list that contains smoothed DO concs
				t.sDO <- t.sDO[,!names(t.sDO)%in%c("doy","year")] # remove the columns that are the year/ doy
				attr(answer, "smoothDO.vec") <- c(t(t.sDO)) # provide the smoothed DO as a long vector, instead of each row containing nobs+2 columns of smoothed DO from a given day
			}
			attr(answer, a0.names[i]) <- answer0[[a0.names[i]]] # assign non "metab" list element as attribute
		}
		
	}else{ # if the list only has one element, or if the list is really a data.frame, then answer0 is what we want
		answer <- answer0
	}
	
	return(answer)
}


