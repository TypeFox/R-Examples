require(RgoogleMaps)

markerdata <- read.csv(system.file('extdata', 'Napoleon.csv', package = 'MSG'))
boundingbox <- qbbox(lon = markerdata[, 'lon'], lat = markerdata[, 'lat'])
NPMap <- GetMap.bbox(boundingbox$lonR, boundingbox$latR, destfile = tempfile(),
                     maptype = 'map', zoom = 6, size = c(640, 370))

lwd_m <- rep(1e-04, nrow(markerdata) - 1)
lwd_m[1:5] <- (15:11) * 1e-5

color <- c('#8B000078', '#00008B78')

tmp <- c(16, 35, 41, 43, 46)
n <- 1 # Engage 1 and Lose 2

PlotOnStaticMap(
  NPMap, FUN = lines,
  lon = markerdata$lon[1:2], lat = markerdata$lat[1:2], 
  lwd = lwd_m[1] * markerdata$size[1], col = color[1], verbose = 0
)

for (i in 2:(nrow(markerdata) - 1)) {
  if (!(i %in% tmp)) {
    PlotOnStaticMap(
      NPMap, FUN = lines,
      lon = markerdata$lon[i:(i + 1)], lat = markerdata$lat[i:(i + 1)],
      lwd = lwd_m[i] * markerdata$size[i], col = color[n], add = TRUE, verbose = 0
    )
  } else {
    n <- ifelse(n == 1, 2, 1)
  }
} 
