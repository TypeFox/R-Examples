
require(SoDA)

load("Examples/with_gps1.rda")

pdf("Examples/gpsScatterPlot.pdf", width = 4, height=4)

with(gps1, {
    par(mar = c(2,2,0,0))
    xy = geoXY(latitude, longitude)
    x = xy[,1]; y = xy[,2]
    plot(x, y, asp = 1)
##   browser()

})

dev.off()
