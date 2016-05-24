cat("\ntest spaMM.filled.contour:")
# spaMM.filled.contour

spaMM.filled.contour(volcano, color = spaMM.colors) # simple

## Comparing the layout with that of filled.contour:
#  (except that it does not achieve the intended effect 
#  in RStudio Plots panel). 

x <- 10*1:nrow(volcano)
y <- 10*1:ncol(volcano)
spaMM.filled.contour(x, y, volcano, color = terrain.colors,
                     plot.title = title(main = "The Topography of Maunga Whau",
                                        xlab = "Meters North", ylab = "Meters West"),
                     plot.axes = { axis(1, seq(100, 800, by = 100))
                                   axis(2, seq(100, 600, by = 100)) },
                     key.title = title(main = "Height\n(meters)"),
                     key.axes = axis(4, seq(90, 190, by = 10)))  # maybe also asp = 1
mtext(paste("spaMM.filled.contour(.) from", R.version.string),
      side = 1, line = 4, adj = 1, cex = .66)

## compare with      

filled.contour(x, y, volcano, color = terrain.colors,
               plot.title = title(main = "The Topography of Maunga Whau",
                                  xlab = "Meters North", ylab = "Meters West"),
               plot.axes = { axis(1, seq(100, 800, by = 100))
                             axis(2, seq(100, 600, by = 100)) },
               key.title = title(main = "Height\n(meters)"),
               key.axes = axis(4, seq(90, 190, by = 10)))  # maybe also asp = 1
mtext(paste("filled.contour(.) from", R.version.string),
      side = 1, line = 4, adj = 1, cex = .66)

# syntax check, no expect_ yet