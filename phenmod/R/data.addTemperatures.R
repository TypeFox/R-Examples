data.addTemperatures <- function(dataset, grid.related.to.Temperatures,
					temperature.filenames, temperature.matrix, 
					temperature.scale.factor,
					out2File=FALSE, silent=FALSE){

	if (is.null(temperature.filenames)) { 
		useTemperatureMatrix <- TRUE 
	} else { 
		useTemperatureMatrix <- FALSE 
	}

	searchFile <- function(files, year){
		for (file in files){
			test <- which(!is.na(match(unlist(strsplit(file,"")),
				unlist(strsplit(as.character(year),"")))))
			check <- seq(from=test[1], by=1, length=nchar(year))
			comparison <- test==check
			if (length(which(comparison==FALSE))==0){
				return(file)
			}
		}
		return("")
	}

	all.years <- sort(unique(dataset$OBS_YEAR_LC))
	edk.one.year <- NA
	# create matrix for temperatures
	temperature <- matrix(NA, length(dataset$STAT_ID), 365)

	# process all years
	temperatures.loaded <- c()
	# prepare Output-Message
	msg <- ""
	for (year in all.years){
		new.year <- TRUE
		if (out2File){
			cat("Year ",year," .. ",sep="")
		}
		observations <- which(dataset$OBS_YEAR_LC == year)
		# iterate over all observations of the processed year
		for (observation.nr in observations){
			# load Temperatures
			if (length(which(temperatures.loaded==year))==0){
				if (useTemperatureMatrix){
					temperature.previous.year <- temperature.matrix[as.character(year),,]
				} else {
					filename <- searchFile(temperature.filenames, year)
					if (file.exists(filename)){
						load(file)
						temperature.previous.year <- edk.one.year
					} else {
						temperature.previous.year <- matrix(NA, 366, 
								dim(grid.related.to.Temperatures)[1])
						warning("Temperature file missing.")
					}
				}

	
				temperatures.loaded <- c(temperatures.loaded, year)
			} else {
				if (new.year){
					new.year <- FALSE
					temperature.previous.year <- temperature.year
				}
			}
			if (length(which(temperatures.loaded==(year+1)))==0){
				if (useTemperatureMatrix){
					temperature.year <- temperature.matrix[as.character(year+1),,]
				} else {
					filename <- searchFile(temperature.filenames, year+1)
					if (file.exists(filename)){
						load(file)
						temperature.year <- edk.one.year
					} else {
						temperature.year <- matrix(NA, 366, 
								dim(grid.related.to.Temperatures)[1])
						warning("Temperature file missing.")
					}
				}

				temperatures.loaded <- c(temperatures.loaded, year+1)
			}

			# add Temperatures to Matrix
			temperature[observation.nr,] <- data.loadTemperature(
						year=dataset$OBS_YEAR_BB[observation.nr],
						temperature.year=temperature.year, 
						temperature.previous.year=temperature.previous.year,
						from.previous.year.doy=dataset$OBS_DAY_LC[observation.nr], 
						length=365,
						position=data.coordinates2gridcellnumber(
							grid=grid.related.to.Temperatures,
							x=dataset$gk4x[observation.nr],
							y=dataset$gk4y[observation.nr]),
						scale.factor=temperature.scale.factor)
			if (!silent){
				if (!out2File){
					cat(rep("\b",nchar(msg)),sep="")
				} else {
					cat("\n",sep="")
				}
				msg <- paste("Year ",year,": ",which(observation.nr == observations),
							" of ", length(observations)," done!",sep="")
				cat(msg,sep="")
				flush.console()
			}
		}
		if (!silent){
			if (out2File){
				cat("Done!\n")
			} else {
				cat("\n")
			}
			msg <- ""
		}
	}

	dataset.with.temperatures <- data.frame(station.id=dataset$STAT_ID, 
							longitude=dataset$STAT_LON, 
							latitude=dataset$STAT_LAT, 
							gk4.x=dataset$gk4x, 
							gk4.y=dataset$gk4y, 
							year.lc=dataset$OBS_YEAR_LC, 
							doy.lc=dataset$OBS_DAY_LC,
							outlier.lc=dataset$OUTLIER_LC, 
							year.bb=dataset$OBS_YEAR_BB, 
							doy.bb.observed=dataset$OBS_DAY_BB, 
							outlier.bb=dataset$OUTLIER_BB,
							temperature=temperature)

	return(dataset.with.temperatures)
}
