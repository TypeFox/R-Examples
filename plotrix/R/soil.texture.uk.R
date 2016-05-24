# UK soil texture plot
soil.texture.uk <- function (soiltexture = NULL, main = "",
 at = seq(0.1, 0.9, by = 0.1),
 axis.labels = c("percent sand", "percent silt", "percent clay"),
 tick.labels = list(l = seq(10, 90, by = 10), r = seq(10, 90, by = 10),
 b = seq(10, 90, by = 10)), show.names = TRUE,
 show.lines = TRUE, col.names = "gray", bg.names = par("bg"),
 show.grid = FALSE, col.axis = "black", col.lines = "gray",
 col.grid = "gray", lty.grid = 3, show.legend = FALSE, label.points = FALSE,
 point.labels = NULL, col.symbols = "black", pch = par("pch"),
 h1 = NA, h3 = NA, t1 = NA, t3 = NA, lwduk = 2, xpos = NA, ypos = NA,
 snames = NA, cexuk = 1.1, ...) {
 
 if(is.na(h1[1])) h1<-c(82, 85, 70, 50, 45, 20) / 100
 if(is.na(h3[1])) h3<-c(18, 15, 30, 30, 35, 0) / 100
 if(is.na(t1[1])) t1<-c(0, 70, 50, 45, 0, 20) / 100
 if(is.na(t3[1])) t3<-c(18, 0, 30, 35, 35, 35) / 100
 # Name positions (x and y, x starting form left point)
 if(is.na(xpos[1])) xpos<-c(0.5,0.77,0.45,0.1,0.45,0.85)
 if(is.na(ypos[1])) ypos<-c(0.65,0.265,0.265,0.07,0.1,0.1)
 if(is.na(snames[1])) snames <- c("Clays","Medium silts","Medium loams",
  "Sands","Light loams","Light silts")
 par(xpd = TRUE)
 plot(0.5, type = "n", axes = FALSE, xlim = c(0,1),ylim = c(0,1),
  main = NA, xlab = NA, ylab = NA)
 triax.plot(x=NULL,main = main, at = at, axis.labels = axis.labels,
  tick.labels = tick.labels, col.axis = col.axis, show.grid = show.grid,
   col.grid = col.grid, lty.grid = lty.grid)
 arrows(0.12, 0.41, 0.22, 0.57, length = 0.15)
 arrows(0.78, 0.57, 0.88, 0.41, length = 0.15)
 arrows(0.6, -0.1, 0.38, -0.1, length = 0.15)
 if(show.lines) {
  triax.segments <- function(h1, h3, t1, t3, col, lwd) {
   segments(1 - h1 - h3/2, h3 * sin(pi/3), 1 - t1 -
   t3/2, t3 * sin(pi/3), col = col, lwd = lwd)
  }
  triax.segments(h1 , h3, t1, t3, col.lines, lwduk)
 }
 if (show.names) {
  boxed.labels(xpos, ypos* sin(pi/3), snames, border = FALSE,
  xpad = 0.5, cex = cexuk)
 }
 par(xpd = FALSE)
 if (is.null(soiltexture)) return(NULL)
 soilpoints <- triax.points(soiltexture, show.legend = show.legend,
  label.points = label.points, point.labels = point.labels,
  col.symbols = col.symbols, pch = pch, ...)
 invisible(soilpoints)
}
