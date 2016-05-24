#' @title Map DAS
#' @description Map sightings, effort, and beaufort from DAS file
#'
#' @param x filename of a DAS file or a \code{data.frame} from \code{\link{das.read}}.
#' @param main main title for plot.
#' @param spp vector of species codes to plot sightings of. If \code{NULL}, only effort will be shown.
#' @param spp.col vector of colors for each species.
#' @param lat.range,lon.range latitude and longitude range of map.
#' @param n.ticks number of desired tick marks for latitude and longitude.
#' @param spp.legend.loc,effort.legend.loc the location of the species and effort legends. 
#'  Can be a single keyword as specified in \code{\link{legend}}, or a two-element numeric vector 
#'  giving the x (longitude) and y (latitude) coordinates.
#'
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @importFrom grDevices rainbow
#' @importFrom graphics abline title lines legend points
#' @importFrom maps map
#' @import mapdata
#' @export
#' 
das.map <- function(x, main, spp = NULL,
                    spp.col = rainbow(length(spp), end = max(1, length(spp) - 2) / length(spp)),
                    lat.range = NULL, lon.range = NULL, n.ticks = 5,
                    spp.legend.loc = "topleft", effort.legend.loc = "topright") {
  das.df <- if(is.data.frame(x)) x else das.read(x)
  das.df <- das.df[!is.na(das.df$Lat) & !is.na(das.df$Long), ]

  # format coordinates
  if(is.null(lat.range)) lat.range <- range(pretty(das.df$Lat))
  if(is.null(lon.range)) lon.range <- range(pretty(das.df$Long))
  lon.range <- ifelse(lon.range < 0, 360 + lon.range, lon.range)
  das.df$Long <- ifelse(das.df$Long < 0, 360 + das.df$Long, das.df$Long)
  in.lat <- das.df$Lat >= min(lat.range) & das.df$Lat <= max(lat.range)
  in.long <- das.df$Long >= min(lon.range) & das.df$Long <= max(lon.range)
  das.df <- das.df[in.lat & in.long, ]

  bft.col <- c("steelblue", "blue", "lightblue", "green", "orange", "red", "darkred")

  # create base map
  map("world2Hires", xlim = lon.range, ylim = lat.range, fill = FALSE)
  abline(h = 0, col = "lightsteelblue")
  title(main = main, line = 3)
  lat.lon.axes(n.ticks)

  # plot effort
  for(i in 1:(nrow(das.df) - 1)) {
    if(das.df$OnEffort[i]) {
      lines(das.df[c(i, i + 1), c("Long", "Lat")],
            lwd = 2,
            col = bft.col[das.df[i, "Bft"] + 1])
    }
  }
  # display effort legend
  effort.x <- effort.legend.loc[1]
  effort.y <- if(is.numeric(effort.legend.loc)) {
    ifelse(effort.legend.loc[2] < 0, 360 + effort.legend.loc[2], effort.legend.loc[2])
  } else NULL
  legend(x = effort.x, y = effort.y, legend = 0:6, lty = "solid", col = bft.col, title = "Beaufort")

  # plot sightings
  if(!is.null(spp)) {
    sights <- das.df[das.df$Event == "A", ]
    spp <- spp[spp %in% unlist(sights[, c("Spp1", "Spp2", "Spp3")])]
    if(length(spp) > 0) {
      for (i in 1:length(spp)) {
        spp.df <- sights[sights$Spp1 == spp[i] | sights$Spp2 == spp[i] | sights$Spp3 == spp[i], ]
        points(spp.df$Long, spp.df$Lat, cex = 1.5, pch = 8, col = spp.col[i])
      }
      # display species legend
      spp.x <- spp.legend.loc[1]
      spp.y <- if(is.numeric(spp.legend.loc)) {
        ifelse(spp.legend.loc[2] < 0, 360 + spp.legend.loc[2], spp.legend.loc[2])
      } else NULL
      legend(x = spp.x, y = spp.y, spp, pch = 8, col = spp.col, horiz = FALSE)
    }
  }

  invisible(das.df)
}
