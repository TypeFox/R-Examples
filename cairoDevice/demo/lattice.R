library(cairoDevice)
Cairo()

require(lattice)
require(stats)

Depth <- equal.count(quakes$depth, number = 8, overlap = 0.1)
xyplot(lat ~ long | Depth, data = quakes)
update(trellis.last.object(), aspect = "iso")

EE <- equal.count(ethanol$E, number = 9, overlap = 1/4)
xyplot(NOx ~ C | EE, data = ethanol, prepanel = function(x, y) prepanel.loess(x, y, span = 1),
            xlab = "Compression Ratio", ylab = "NOx (micrograms/J)",
            panel = function(x, y) { panel.grid(h=-1, v= 2); panel.xyplot(x, y); panel.loess(x,y, span=1); },
            aspect = "xy")
			
plot <- xyplot(sunspot.year ~ 1700:1988, xlab = "", type = "l",
                    scales = list(x = list(alternating = 2)),
                    main = "Yearly Sunspots")
     print(plot, position = c(0, .3, 1, .9), more = TRUE)
     print(update(plot, aspect = "xy", main = "", xlab = "Year"),
           position = c(0, 0, 1, .3))


