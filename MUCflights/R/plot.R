
# wie kann plot.routes die Argumente von movie.routes automatisch uebernehmen?
plot.routes <- function(x, time,
                        borders, capitals,
                        col_flights = c(S = "#2F984F", L = "#2E7EBC"),
                        start = c(lon = 11.78609, lat = 48.35378),
                        citynames = FALSE, NightDay = FALSE, ...){

  drawMap(borders, capitals, ...)

  if (NightDay) plot(NightDay(time, timezone = 1), add = TRUE)

  points(start["lon"], start["lat"], col = 2, pch = 19)
  mtext(as.character(time), 3, line = 0)

         a <- function(i) {
          if(any(tmploc$lat[i] > (start["lat"] - 1) & tmploc$lat[i] < (start["lat"] + 1) &
             tmploc$lon[i] > (start["lon"] - 1) & tmploc$lon[i] < (start["lon"] + 1))) {
              points(tmploc$Longitude[i], tmploc$Latitude[i], col = "red", pch = 19, cex = 0.5)
              if(citynames == TRUE)
                  text(tmploc$Longitude[i], tmploc$Latitude[i], tmploc$cities[i])
          }
          ax <- gcIntermediate(tmploc[i, c("lon", "lat")],
                               tmploc[i, c("Longitude", "Latitude")], addStartEnd = TRUE)
          ay <- matrix(ncol = 2, c(tmploc[i, c("lon", "lat")]), byrow = TRUE)
          lsk.i <- unique(tmploc[i, "lsk"])
          cl <- col_flights[lsk.i]
          points(ay, col = cl, pch = 19, cex = 0.5)
          lines(ax, col = "grey70", lty = 3)
      }

  #tmploc <- subset(x, x$time == time)
  tmploc <- x[x$time == time, ]
  if ( nrow(tmploc) >= 1 ) {
    tapply(1:nrow(tmploc), tmploc$fnr[, drop = TRUE], a)
  }
}

drawMap <- function(borders, capitals, col_border = "darkgrey",
                    col_land = "#FFFFE5", bg = "#F7FBFF", axes = TRUE, ...) {

  plot(borders, border = col_border, col = col_land, bg = bg, axes = axes)
  box()

  if(!missing(capitals)){
    points(capitals, ...)
  }
}


###### plotting examples
#system.time(drawMap(shapeEU, spCaps, bg="#F7FBFF"))
#drawMap(shapeEU, spCaps, bg="#F7FBFF", border="darkgrey")
## main, xlab, ylab etc. scheinen nicht zu funktionieren
