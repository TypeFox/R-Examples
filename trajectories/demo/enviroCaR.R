## Shows how to se the enviroCaR R package to import enviroCar tracks.

if (!require(devtools))
	stop("install devtools first")

if (!require(enviroCaR)) {
	devtools::install_github("enviroCar/enviroCaR"); 
	library(enviroCaR)
}

## get all track ids
ids <- getTrackIDs("https://envirocar.org/api/stable")

## get ids using a bounding box and a time interval
## bounding box
bbox <- matrix(c(7.318136,51.802163, 7.928939,52.105665), nrow=2, ncol=2)

## time interval: 96 hours
t2 <- as.POSIXct("2015-01-20 10:11:20 CEST")
t1 <- t2 - as.difftime(96, unit="hours")

## get ids
ids <- getTrackIDs("https://envirocar.org/api/stable", bbox, list(first.time = t1, last.time = t2))

## import tracks, returns a TracksCollection
trcol <- importEnviroCar("https://envirocar.org/api/stable", ids)

## import single track, returns a Tracks object
track <- importSingleTrack("https://envirocar.org/api/stable", ids[1])
