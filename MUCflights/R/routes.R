
routes <- function(flights.info, start.IATA = "MUC") {
        
        airports <- NULL
        #load(system.file("data/airports.RData", package = "MUCflights"))
        data("airports", package = "MUCflights", envir = sys.frame(which = 1))
        stopifnot(!is.null(airports))

	# Fluege mit Ziel Muenchen loeschen
	flights.info <- subset(flights.info, 
            !((flights.info$lsk == "S") & flights.info$ha1 == "MUC"))

	flights.info <- merge (flights.info, airports, by.x = "ha1", by.y = "IATA")

	# Distanzen von Startflughafen (start) in Kilometer ausrechnen
	R <- 6378137 / 1000
	start <- as.matrix(airports[airports$IATA ==
										start.IATA, c("Longitude", "Latitude")])
	flights.info$distance <- as.numeric(
		distHaversine(flights.info[, c("Longitude", "Latitude")], start, r = R))
		
	flights.info <- subset(flights.info, flights.info$distance > 1)


	# Datensatz sortieren nach stt
	flights.info <- (flights.info[order(flights.info$stt), ])

	# Dauer je Flug in Min. berechen
	duration <- round(flights.info$distance / 900 * 60)

	# Zeiten codieren
	stt <- as.POSIXct(strftime(flights.info$stt, format = "%Y-%m-%d %H:%M:%S.0"))
	ett <- as.POSIXct(strftime(flights.info$ett, format = "%Y-%m-%d %H:%M:%S.0"))

	S.ind <- which(flights.info$lsk == "S")
	L.ind <- which(flights.info$lsk == "L")

	lonlat.S <- flights.info[S.ind, c("Longitude", "Latitude")]
	lonlat.L <- flights.info[L.ind, c("Longitude", "Latitude")]

	## zwischenpunkte der flugroute berechnen
	gcIntermed.S <- do.call(rbind, gcIntermediate(start,
												lonlat.S, n = duration[S.ind]))
	gcIntermed.L <- do.call(rbind, gcIntermediate(lonlat.L,
												start, n = duration[L.ind]))

	## time und delay sequenzen
	## unlist konvertiert die Zeiten leider zu integers,
	## deshalb unten der workaround mit der '<-c("POSIXct", "POSIXt")'-Zuweisung
	time.S <- unlist(lapply(S.ind, function(z)
				seq.POSIXt(from = stt[z], by = 60, length.out = duration[z])))

	time.L <- unlist(lapply(L.ind, function(z)
			rev(seq.POSIXt(from = stt[z], by = -60, length.out = duration[z]))))

	delay.S <- unlist(lapply(S.ind, function(z)
				seq.POSIXt(from = ett[z], by = 60, length.out = duration[z])))

	delay.L <- unlist(lapply(L.ind, function(z)
			rev(seq.POSIXt(from = ett[z], by = -60, length.out = duration[z]))))

	timeANDdelay <- rbind(data.frame(time = time.S, delay = delay.S),
									data.frame(time = time.L, delay = delay.L))

	for (i in 1:ncol(timeANDdelay))
            class(timeANDdelay[, i]) <- c("POSIXct", "POSIXt")

	## set id
	id.S <- rep((1:nrow(flights.info))[S.ind], times = duration[S.ind])
	id.L <- rep((1:nrow(flights.info))[L.ind], times = duration[L.ind])

	## datensatz zusammenbauen
	posANDtime <- cbind(id = c(id.S, id.L),
								rbind(gcIntermed.S, gcIntermed.L), timeANDdelay)

	## prepare flights.info
	flights.info$stt <- stt
	flights.info$ett <- ett
	flights.info$duration <- duration
	flights.info$id <- 1:nrow(flights.info)

	## merge position-matrix and flight.info
	positionMatrix <- merge(posANDtime, flights.info, by = "id")
	positionMatrix[order(positionMatrix$id),]
	# Klasse festlegen
	class(positionMatrix) <- c("routes", class(positionMatrix))
	return(positionMatrix)

}
