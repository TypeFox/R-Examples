data.combine.stationNet <- function(dataset, range, 
					alt.range, silent=FALSE, 
					out2File=FALSE){
	getDistance <- function(x.1, y.1, x.2, y.2){
		distance <- sqrt((x.1-x.2)^2 + (y.1-y.2)^2)
		return(distance)
	}
	getCoordinates <- function(station, dataset){
		station.nr <- which(dataset$STAT_ID==station)[1]
		return(c(dataset$gk4x[station.nr],dataset$gk4y[station.nr]))
	} 
	getAltitude <- function(station, dataset){
		station.nr <- which(dataset$STAT_ID==station)[1]
		return(dataset$STAT_ALT[station.nr])
	} 
	changeOrder <- function(vec, pos1, pos2){
		tmp <- vec[pos1]
		vec[pos1] <- vec[pos2]
		vec[pos2] <- tmp
		return(vec)
	}

	stations <- as.list(unique(dataset$STAT_ID))
	stations.net <- list()

	## output message
	msg <- ""

	## iterate over all stations
	for (station in stations){
		## get station coordinates
		station.coordinates <- getCoordinates(station, dataset)
		station.altitude <- getAltitude(station, dataset)

		## search near stations
		near.stations <- c()
		distances <- c()
		
		## iterate over all other stations
		for (station.to.check in stations){
			if (station.to.check != station){
				## get station coordinates
				coordinates.to.check <- getCoordinates(station.to.check, dataset)
				## get distance to station of main iteration
				distance <- getDistance(station.coordinates[1],station.coordinates[2],
						coordinates.to.check[1], coordinates.to.check[2])
				## get altitude difference
				altitude <- getAltitude(station.to.check, dataset)
				altitude.dif <- abs(altitude-station.altitude)
				rm(altitude)
				## check if distance is lower than range
				## and if altitude is in alt.range
				if ((distance <= range)&&(altitude.dif <= alt.range)){
					## if true, set station as a near station
					near.stations[length(near.stations)+1] <- station.to.check
					distances[length(distances)+1] <- distance
				}
			}
		}
		if (length(distances) > 0){
			## order stations, station with lowest distance as first
			distances.sorted <- sort(distances)
			near.stations.sorted <- c()
	
			count <- 0
			max.count <- 0
			already.set <- FALSE		

			for (distance in distances.sorted){
				length.same.distances <- length(which(distance == distances))
				if (length.same.distances == 1){
					distance.id <- which(distance == distances)
				} else {
					if (!already.set){
						max.count <- length(which(distance == distances))
						already.set <- TRUE
						count <- 1
					}
					distance.id <- which(distance == distances)[count]
			
					count <- count + 1
					if (count > max.count){
						already.set <- FALSE
					}
				}

				near.stations.sorted[length(near.stations.sorted)+1] <- 	
					near.stations[distance.id]

			}
			near.stations <- near.stations.sorted
			stations.net[[length(stations.net)+1]] <- near.stations
			rm(distances, distances.sorted, near.stations.sorted, near.stations)
		} else {
			stations.net[[length(stations.net)+1]] <- NA
			rm(distances, near.stations)
		}
		## output message
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
		cat("\n",sep="")
	}
	return(stations.net)
}