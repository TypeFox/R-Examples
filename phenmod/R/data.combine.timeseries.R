data.combine.timeseries <- function(dataset, clusters, 
					silent=FALSE, out2File=FALSE, 
					minimalClusterSize=5){

	# vectors which should be combined to data.frame
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

	# output msg
	msg <- ""
	count <- 1

	# iterate over all clusters
	for (cluster in clusters){
		if (length(cluster) >= minimalClusterSize){
			success <- TRUE
			stations <- cluster
			station.numbers <- c()
			for (station in stations){
				station.numbers <- c(station.numbers, 
					which(dataset$STAT_ID == station))
			}
			data.tmp <- dataset[station.numbers,]
			## remove outliers
			data.tmp <- data.tmp[which(data.tmp$OUTLIER_LC==0),]
			data.tmp <- data.tmp[which(data.tmp$OUTLIER_BB==0),]
			# extract values
			x.bb <- data.tmp$OBS_DAY_BB
			x.lc <- data.tmp$OBS_DAY_LC
			f1.bb <- data.tmp$OBS_YEAR_BB
			f1.lc <- data.tmp$OBS_YEAR_LC
			f2 <- data.tmp$STAT_ID
		
			data.to.combine.bb <- data.frame(x=x.bb, f1=f1.bb, f2=f2)
			pheno.fit.bb <- try(pheno.lad.fit(data.to.combine.bb), silent=TRUE)
			if (!inherits(pheno.fit.bb, "try-error")){
				obsday.bb <- pheno.fit.bb$f1
				obsyear.bb <- as.integer(pheno.fit.bb$f1.lev)
			} else {
				success <- FALSE
			}

			data.to.combine.lc <- data.frame(x=x.lc, f1=f1.lc, f2=f2)
			pheno.fit.lc <- try(pheno.lad.fit(data.to.combine.lc), silent=TRUE)

			if (!inherits(pheno.fit.lc, "try-error")){
				obsday.lc <- pheno.fit.lc$f1
				obsyear.lc <- as.integer(pheno.fit.lc$f1.lev)
			} else {
				success <- FALSE
			}

			if (success){
				# coordinates mean
				lon.values <- c()
				lat.values <- c()
				gk4x.values <- c()
				gk4y.values <- c()

				for (station in stations){
					lon.values[length(lon.values)+1] <- 
						(data.tmp[which(data.tmp$STAT_ID == station),])$STAT_LON[1]
					lat.values[length(lat.values)+1] <- 
						(data.tmp[which(data.tmp$STAT_ID == station),])$STAT_LAT[1]
					gk4x.values[length(gk4x.values)+1] <- 
						(data.tmp[which(data.tmp$STAT_ID == station),])$gk4x[1]
					gk4y.values[length(gk4y.values)+1] <- 
						(data.tmp[which(data.tmp$STAT_ID == station),])$gk4y[1]
				
				}

				stat.lon <- mean(lon.values, na.rm=TRUE)
				stat.lat <- mean(lat.values, na.rm=TRUE)
				gk4x.value <- mean(gk4x.values, na.rm=TRUE)
				gk4y.value <- mean(gk4y.values, na.rm=TRUE)

				# alt mean
				alt.values <- data.tmp$STAT_ALT
				stat.alt <- mean(alt.values, na.rm=TRUE)

				rm(lon.values, lat.values, gk4x.values, gk4y.values, alt.values)

				# new station id
				stat.id <- ""
				for (station in stations){
					if (nchar(stat.id)==0){
						stat.id <- paste(station,sep="")
					} else {
						stat.id <- paste(stat.id, station, sep=",")
					}
				}
	
				# add to vectors
				STAT_ID <- c(STAT_ID, rep(stat.id, length(obsday.bb)))
				STAT_LON <- c(STAT_LON, rep(stat.lon, length(obsday.bb)))
				STAT_LAT <- c(STAT_LAT, rep(stat.lat, length(obsday.bb)))
				STAT_ALT <- c(STAT_ALT, rep(stat.alt, length(obsday.bb)))
				OBS_YEAR_LC <- c(OBS_YEAR_LC, obsyear.lc)
				OBS_DAY_LC <- c(OBS_DAY_LC, obsday.lc)
				OUTLIER_LC <- c(OUTLIER_LC, rep(0, length(obsday.bb)))
				OBS_YEAR_BB <- c(OBS_YEAR_BB, obsyear.bb)
				OBS_DAY_BB <- c(OBS_DAY_BB, obsday.bb)
				OUTLIER_BB <- c(OUTLIER_BB, rep(0, length(obsday.bb)))
				gk4x <- c(gk4x, rep(gk4x.value, length(obsday.bb)))
				gk4y <- c(gk4y, rep(gk4y.value, length(obsday.bb)))
			}
		}
		if (!silent){
			if (out2File){
				cat("\n")
			} else {
				cat(rep("\b",nchar(msg)),sep="")
			}			
			msg <- paste(count, " of ", length(clusters), " clusters done!",sep="")
			cat(msg)
			flush.console()
			count <- count+1
		}
	}

	if (!silent){
		cat("\n",sep="")
	}
	
	# create dataframe
	dataset.combine <- data.frame(STAT_ID=STAT_ID, STAT_LON=STAT_LON, 
				STAT_LAT=STAT_LAT, STAT_ALT=STAT_ALT, 
				OBS_YEAR_LC=OBS_YEAR_LC, OBS_DAY_LC=OBS_DAY_LC,
				OUTLIER_LC=OUTLIER_LC, OBS_YEAR_BB=OBS_YEAR_BB, 
				OBS_DAY_BB=OBS_DAY_BB, OUTLIER_BB=OUTLIER_BB, 
				gk4x=gk4x, gk4y=gk4y)

	return(dataset.combine)
}