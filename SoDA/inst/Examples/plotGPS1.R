source("Examples/trackArrows.R")

source("Examples/plotGPSArrows.R")

require(SoDA)

load("Examples/with_gps1.rda")

pdf("Examples/plotGPS1.pdf", width = 4, height=4)

object = new("GPSTrack", latitude = gps1$latitude,
  longitude = gps1$longitude, elevation = gps1$elevation,
  time = gps1$date)

    par(mar = c(2,2,1,.1))
    plotGPSArrows(object)
##   browser()



dev.off()
