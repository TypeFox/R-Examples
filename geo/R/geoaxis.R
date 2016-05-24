#' Setup axes for a geoplot
#' 
#' Set up axes for a geoplot.
#' 
#' 
#' @param side On which side should the axis be drawn?
#' @param pos Positions of labels (?)
#' @param dist Distance from plot of axes labels (?)
#' @param dlat Latitude resolution of labels (?)
#' @param dlon Longitude resolution of labels (?)
#' @param csi Character size inches, default 0.12
#' @param cex Character size expansion, default 0.7
#' @param inside Whether or not stay inside ??????????
#' @param r yet another setting ?????????
#' @param \dots Additional arguments to \code{geotext}
#' @return Adds labels to a geoplot, no value returned.
#' @note Needs further elaboration.
#' @seealso Calls \code{\link{geotext}}, called by \code{\link{init}}.
#' @keywords aplot
#' @export geoaxis
geoaxis <-
function(side, pos, dist, dlat = 0.5, dlon = 1,csi=0.12, cex = 0.7, inside = T, r = 1,
	...)
{
	geopar <- getOption("geopar")
	m <- par()$cex * csi
	if(inside)
		m <-  - m
	if(side == 2 || side == 4)
		ratio <- diff(geopar$origin$lon)/geopar$gpar$pin[1]
	else ratio <- diff(geopar$origin$lat)/geopar$gpar$pin[2]
	if(side == 2 || side == 4) {
		if(missing(pos)) {
			pos1 <- (geopar$origin$lat[1] %/% dlat) * dlat - dlat
			pos2 <- (geopar$origin$lat[2] %/% dlat) * dlat + dlat
			pos <- seq(pos1, pos2, by = dlat)
		}
		pos <- pos[pos <= geopar$origin$lat[2] & pos >= geopar$origin$
			lat[1]]
		if(missing(dist)) {
			if((side == 4 && inside) || (side == 2 && !inside)) {
				lat1 <- pos %% 1
				lat2 <- lat1 * 60 %% 1
				if(any(lat2))
					lm <- 8
				else if(any(lat1))
					lm <- 6
				else lm <- 4
				dist <- lm * m * 0.6 * r
			}
			else if((side == 2 && inside) || (side == 4 && !inside)
				)
				dist <- m/2.5 * r
		}
	}
	else if(side == 3 || side == 1) {
		if(missing(dist))
			dist <- m * r
		if(missing(pos)) {
			pos1 <- (geopar$origin$lon[1] %/% dlon) * dlon - dlon
			pos2 <- (geopar$origin$lon[2] %/% dlon) * dlon + dlon
			pos <- seq(pos1, pos2, by = dlon)
		}
		pos <- pos[pos <= geopar$origin$lon[2] & pos >= geopar$origin$
			lon[1]]
	}
	if(side == 2 || side == 4) {
		if(side == 2)
			lonpos <- geopar$origin$lon[1] - ratio * dist
		if(side == 4)
			lonpos <- geopar$origin$lon[2] + ratio * dist
		lat1 <- trunc(pos)
		lat2 <- pos %% 1
		txt <- paste(lat1, "\u00b0", sep = "")
		lat2 <- round(lat2 * 60, 2)
		i <- lat2 > 0
		if(any(i))
			txt[i] <- paste(txt[i], lat2[i], "'", sep = "")
		geotext(pos, rep(lonpos, length(lat1)), adj = 0, txt, cex = cex,
			outside = T, ...)
	}
	else {
		if(side == 1)
			latpos <- geopar$origin$lat[1] - ratio * dist
		if(side == 3)
			latpos <- geopar$origin$lat[2] + ratio * dist
		pos1 <- abs(pos)
		lon1 <- trunc(pos1)
		lon2 <- pos1 %% 1
		txt <- paste(lon1, "\u00b0", sep = "")
		lon2 <- round(lon2 * 60, 2)
		i <- lon2 > 0
		if(any(i))
			txt[i] <- paste(txt[i], lon2[i], "'", sep = "")
		geotext(rep(latpos, length(pos)), pos, txt, adj = 0.5, cex = 
			cex, outside = T, ...)
	}
}

