data.extract <- function(data.budburst, data.leafcolouring, 
			valid.years=1952:2009, out2File=FALSE, 
			silent=FALSE){
	# extract stations
	stations.bb <- unique(data.budburst$DWD_STAT_ID)
	stations.lc <- unique(data.leafcolouring$DWD_STAT_ID)

	# create combined dataset vectors
	## all observing stations
	stations <- unique(c(stations.bb, stations.lc))
	
	## vectors which should be combined to data.frame
	STAT_ID <- c()
	STAT_LON <- c()
	STAT_LAT <- c()
	STAT_ALT <- c()
	OBS_YEAR_BB <- c()
	OBS_DAY_BB <- c()
	OUTLIER_BB <- c()
	OBS_YEAR_LC <- c()
	OBS_DAY_LC <- c()
	OUTLIER_LC <- c()
	gk4x <- c()
	gk4y <- c()

	# create output message
	msg <- ""

	# iterate over all stations
	for (station in stations){
		# check if station observes BB
		if (length(which(stations.bb == station)) > 0){
			station.observes.bb <- TRUE
		} else {
			station.observes.bb <- FALSE
		}
		# check if station observes LC
		if (length(which(stations.lc == station)) > 0){
			station.observes.lc <- TRUE
		} else {
			station.observes.lc <- FALSE
		}

		# only use stations which observed a lc-doy
		if (station.observes.lc){
			selected <- which(data.leafcolouring$DWD_STAT_ID == station)
			# iterate over all datapoints of that station
			for (select in selected){
				year.lc <- data.leafcolouring[select,,]$OBS_YEAR
				station_longitude <- data.leafcolouring[select,,]$STAT_LON
				station_latitude <- data.leafcolouring[select,,]$STAT_LAT
				station_altitude <- data.leafcolouring[select,,]$STAT_ALT
				station_name <- data.leafcolouring[select,,]$STAT_NAME
				doy.lc <- data.leafcolouring[select,,]$OBS_DAY
				outlier.lc <- data.leafcolouring[select,,]$outlier
				# search for budbust-doy
				year.bb <- year.lc + 1
				doy.bb <- NA
				outlier.bb <- NA
				if (station.observes.bb){
					related.bb.stations.nr <- 
						which(data.budburst$DWD_STAT_ID == station)
					related.bb.stations <- 
						data.budburst[related.bb.stations.nr,,]
					dataset.nr.with.same.year <- 
						which(related.bb.stations$OBS_YEAR == year.bb)
					datasets.with.same.year <- 
						related.bb.stations[dataset.nr.with.same.year,,]
					doy.bb <- datasets.with.same.year$OBS_DAY
					outlier.bb <- datasets.with.same.year$outlier
					if (length(doy.bb) == 0){
						doy.bb <- NA
						outlier.bb <- NA
					}
					if (length(doy.bb) > 1){
						doy.bb <- doy.bb[which(is.na(doy.bb)==FALSE)]
						outlier.bb <- outlier.bb[which(is.na(doy.bb)==FALSE)]
						doy.bb <- doy.bb[1]
						outlier.bb <- outlier.bb[1]
					}
				}
				# add datapoint to vectors
				# but only for valid years
				if (length(which(year.bb == valid.years)) > 0){
					STAT_ID <- c(STAT_ID, station)
					STAT_LON <- c(STAT_LON, station_longitude)
					STAT_LAT <- c(STAT_LAT, station_latitude)
					STAT_ALT <- c(STAT_ALT, station_altitude)
					OBS_YEAR_LC <- c(OBS_YEAR_LC, year.lc)
					OBS_DAY_LC <- c(OBS_DAY_LC, doy.lc)
					OUTLIER_LC <- c(OUTLIER_LC, outlier.lc)
					OBS_YEAR_BB <- c(OBS_YEAR_BB, year.bb)
					OBS_DAY_BB <- c(OBS_DAY_BB, doy.bb)
					OUTLIER_BB <- c(OUTLIER_BB, outlier.bb)
				}
			}
		}

		if (!silent){
			if (out2File){
				cat("\n")
			} else {
				cat(rep("\b",nchar(msg)),sep="")
			}
			msg <- paste(which(stations == station), " of ", 
				length(stations), " Stations done!",sep="")
			cat(msg)
			flush.console()
		}
	}
	if (!silent){
		cat("\n")
	}

	# add Gauss-Krueger-Coordinates
	new.coordinates <- util.geoco2gk(STAT_LON, STAT_LAT, 4)
	gk4x <- new.coordinates[,1]
	gk4y <- new.coordinates[,2]

	dataset <- data.frame(STAT_ID=STAT_ID, STAT_LON=STAT_LON, 
				STAT_LAT=STAT_LAT, STAT_ALT=STAT_ALT, 
				OBS_YEAR_LC=OBS_YEAR_LC, OBS_DAY_LC=OBS_DAY_LC,
				OUTLIER_LC=OUTLIER_LC, OBS_YEAR_BB=OBS_YEAR_BB, 
				OBS_DAY_BB=OBS_DAY_BB, OUTLIER_BB=OUTLIER_BB, 
				gk4x=gk4x, gk4y=gk4y)

	return(dataset)
}