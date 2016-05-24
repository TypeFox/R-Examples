#scalebar


# Calculates the geodesic distance between two points specified by radian latitude/longitude using the
# Haversine formula (hf) from: http://www.r-bloggers.com/great-circle-distance-calculations-in-r/
.torad <- function(deg) {
  deg*pi/180.0
}

.tolatlon <- function(x, y, epsg) {
  rgdal::CRSargs(sp::CRS(paste0("+init=epsg:", epsg))) #hack to make sure rgdal stays in Imports:
  coords <- sp::coordinates(matrix(c(x,y), byrow=TRUE, ncol=2))
  spoints <- sp::SpatialPoints(coords, sp::CRS(paste0("+init=epsg:", epsg)))
  spnew <- sp::spTransform(spoints, sp::CRS("+init=epsg:4326"))
  c(sp::coordinates(spnew)[1], sp::coordinates(spnew)[2])
}

.geodist <- function(x1, y1, x2, y2, epsg) {
  lonlat1 <- .tolatlon(x1, y1, epsg)
  lonlat2 <- .tolatlon(x2, y2, epsg)

  long1 <- .torad(lonlat1[1])
  lat1 <- .torad(lonlat1[2])
  long2 <- .torad(lonlat2[1])
  lat2 <- .torad(lonlat2[2])
  R <- 6371009 # Earth mean radius [m]
  delta.long <- (long2 - long1)
  delta.lat <- (lat2 - lat1)
  a <- sin(delta.lat/2)^2 + cos(lat1) * cos(lat2) * sin(delta.long/2)^2
  c <- 2 * asin(min(1,sqrt(a)))
  d = R * c
  return(d) # Distance in m
}

.fromsi <- function(sivalue, unit) {
  if(unit == "km") {
    sivalue / 1000.0
  } else if(unit == "m") {
    sivalue
  } else if(unit =="ft") {
    sivalue * 3.28084
  } else if(unit == "mi") {
    sivalue / 1609.344051499
  } else if(unit == "in") {
    sivalue * 39.370079999999809672
  } else if(unit == "cm") {
    sivalue * 100.0
  } else {
    stop("Unrecognized unit: ", unit)
  }
}

.tosi <- function(unitvalue, unit) {
  if(unit == "km") {
    unitvalue * 1000.0
  } else if(unit == "m") {
    unitvalue
  } else if(unit =="ft") {
    unitvalue / 3.28084
  } else if(unit == "mi") {
    unitvalue * 1609.344051499
  } else if(unit == "in") {
    unitvalue / 39.370079999999809672
  } else if(unit == "cm") {
    unitvalue / 100.0
  } else {
    stop("Unrecognized unit: ", unit)
  }
}

#' Get Scale Bar Parameters
#'
#' Get default scale bar parameters based on the current plot (i.e. \code{par("usr")}).
#' The algorithm attempts to detect the best equally divisable distance to use for the
#' scale bar, and returns a \code{list} object with attributes that allow any type of
#' scale bar to be drawn. The only way to manipulate the values chosen by the algorithm
#' is to change the \code{widthhint} argument. For generic XY plots, pass \code{plotunit}.
#'
#' @param plotunit The unit which the current plot is plotted in, one of \code{cm},
#' \code{m}, \code{km}, \code{in}, \code{ft}, \code{mi}. or \code{latlon}. This
#' parameter is optional if \code{plotepsg} is passed.
#' @param plotepsg The projection of the current plot. If extents are valid lat/lons,
#' the projection is assumed to be lat/lon (EPSG:4326), or Spherical Mercator otherwise
#' (EPSG:3857). This is done to work seamlessly with OpenStreetMap packages.
#' @param widthhint The fraction of the plottable width which the scale bar should
#' (mostly) occupy.
#' @param unitcategory One of "metric" or "imperial"
#' @return a \code{list} of parameters: \code{$widthu} (width of the scalebar in human
#' readable units); \code{$unit} (the human readable unit); \code{$majordivu} (the size
#' of the divisions in human readable units); \code{$majordivs} (the number of divisions);
#' \code{$widthplotunit} (width of the scalebar in plotting units); \code{$majordivplotunit}
#' (the width of divisions in plotting units); \code{$labeltext} (label text); and \code{extents}
#' the user extents (\code{par('usr')}) that were used to calculate the parameters.
#' @export
#'
#' @examples
#' plot(1:5, 1:5, asp=1)
#' scalebarparams(plotunit="m")
#' \donttest{
#' library(maptools)
#' data(wrld_simpl)
#' plot(wrld_simpl, xlim=c(-66.86, -59.75), ylim=c(43, 47.3)) #Nova Scotia
#' scalebarparams()
#' }
#'
#' @seealso \link{addscalebar}
#'
#'
scalebarparams <- function(plotunit=NULL, plotepsg=NULL, widthhint=0.25, unitcategory="metric") {
  #params check
  if(!(unitcategory %in% c("metric", "imperial"))) stop("Unrecognized unitcategory: ", unitcategory)

  extents <- graphics::par('usr')
  if(is.null(plotepsg) && is.null(plotunit)) {
    #check for valid lat/lon in extents
    if(extents[1] >= -180 &&
       extents[1] <= 180 &&
       extents[2] >= -180 &&
       extents[2] <= 180 &&
       extents[3] >= -90 &&
       extents[3] <= 90 &&
       extents[4] >= -90 &&
       extents[4 <= 90]) {
      warning("Autodetect projection: assuming lat/lon (epsg 4326)")
      plotepsg <- 4326
    } else {
      #else assume google mercator used by {OpenStreetMap} (epsg 3857)
      warning("Audotdetect projection: assuming Google Mercator (epsg 3857)")
      plotepsg <- 3857
    }
  } else if(!is.null(plotunit)) {
    if(plotunit=="latlon") {
      plotepsg <- 4326
    }
  }

  if(!is.null(plotepsg)) {
    widthbottom <- .geodist(extents[1], extents[3], extents[2], extents[3], plotepsg)
    midY <- mean(c(extents[3], extents[4]))
    widthmiddle <- .geodist(extents[1], midY, extents[2], midY, plotepsg)
    widthtop <- .geodist(extents[1], extents[4], extents[2], extents[4], plotepsg)
    percentdiff <- (max(widthbottom, widthmiddle, widthtop) -
                      min(widthbottom, widthmiddle, widthtop)) / min(widthbottom, widthmiddle, widthtop)
    if(percentdiff > .05) warning("Scale on map varies by more than 5%, scalebar may be inaccurate")
    widthm <- widthmiddle
    mperplotunit <- widthmiddle/(extents[2]-extents[1])
  } else {
    heightm <- .tosi(extents[4] - extents[3], plotunit)
    widthm <- .tosi(extents[2] - extents[1], plotunit)
    mperplotunit <- .tosi(1.0, plotunit)
  }

  geowidthm <- widthm * widthhint

  if(geowidthm < 1) {
    scaleunits <- c("cm", "in")
  } else if(geowidthm < 1600) {
    scaleunits <- c("m", "ft")
  } else {
    scaleunits <- c("km", "mi")
  }

#   String unit = units[unitCategory] ;
  if(unitcategory == "metric") {
    unit <- scaleunits[1]
  } else {
    unit <- scaleunits[2]
  }
#   double widthHintU = Units.fromSI(geoWidthM, unit) ;
  widthhintu <- .fromsi(geowidthm, unit)
#   double tenFactor = Math.floor(Math.log10(widthHintU)) ;
  tenfactor <- floor(log10(widthhintu))
#   double widthInTens = Math.floor(widthHintU / Math.pow(10, tenFactor)) ;
  widthintens <- floor(widthhintu / (10^tenfactor))
  if(widthintens == 1) {
    widthintens <- 10
    tenfactor = tenfactor - 1 ;
  } else if(widthintens == 7) {
    widthintens <- 6
  } else if(widthintens == 9) {
    widthintens <- 8
  }

  if(widthintens < 6) {
    majdivtens <- 1
  } else {
    majdivtens <- 2
  }

#   double widthU = widthInTens * Math.pow(10, tenFactor) ;
  widthu <- widthintens * 10^tenfactor
#   double majorDiv = majDivTens * Math.pow(10, tenFactor) ;
  majordiv <- majdivtens * 10^tenfactor
#   long majorDivs = Math.round(widthU / majorDiv) ;
  majordivs <- round(widthu / majordiv)
#   double widthPx = Units.toSI(widthU, unit) / mPerPixel ;
  widthplotunit <- .tosi(widthu, unit) / mperplotunit
#   double majorDivPx = widthPx / majorDivs ;
  majordivplotunit <- widthplotunit / majordivs
#   this.scaleParameters = new double[] {widthU, majorDiv, widthPx, majorDivPx} ;
  params = list()
  params$widthu <- widthu
  params$unit <- unit
  params$majordivu <- majordiv
  params$majordivs <- majordivs
  params$widthplotunit <- widthplotunit
  params$majordivplotunit <- majordivplotunit
  params$labeltext <- paste(as.integer(widthu), unit)
  params$extents <- extents
#   this.labelText = String.valueOf(Math.round(widthU)) + " " + unit ;
  params

}

#' Raw Plot Scale Bar
#'
#' Just in case anybody is hoping to draw a custom scalebar, this is the
#' method used to plot it. If you don't know what this is, you should probably
#' be using \link{addscalebar}.
#'
#' @param x The position (user) to draw the scale bar
#' @param y The position (user) to draw the scale bar
#' @param ht The height(in user coordinates) to draw the scale bar
#' @param params Scalebar parameters as generated by \link{scalebarparams}
#' @param style One of \code{bar} or \code{ticks}
#' @param adj Where to align the scale bar relative to \code{x} and \code{y}
#' @param tick.cex If \code{style=="ticks"}, the height of interior ticks.
#' @param bar.cols A vector of color names to be repeated for a \code{bar}
#' style scalebar.
#' @param lwd Passed when drawing lines associated with the scalebar
#' @param linecol Passed when drawing lines associated with the scalebar
#'
#' @seealso \link{addscalebar}
#'
#' @export
#'
#'
plotscalebar <- function(x, y, ht, params, style="bar", adj=c(0,0), tick.cex=0.7,
                         bar.cols=c("black", "white"), lwd=1, linecol="black") {
  wd <- params$widthplotunit
  if(style=="bar") {
    cols <- rep(bar.cols, params$majordivs/(length(bar.cols))+1)
    for(i in 1:params$majordivs) graphics::rect(x-adj[1]*wd+(i-1)*params$majordivplotunit, y-adj[2]*ht+ht,
                                            x-adj[1]*wd+i*params$majordivplotunit, y-adj[2]*ht, col=cols[i],
                                            lwd=lwd, border=linecol)
  } else if(style=="ticks") {
    outerx <- c(x-adj[1]*wd,
                x-adj[1]*wd,
                x-adj[1]*wd+wd,
                x-adj[1]*wd+wd)
    outery <- c(y-adj[2]*ht+ht,
                y-adj[2]*ht,
                y-adj[2]*ht,
                y-adj[2]*ht+ht)
    graphics::lines(outerx, outery, lwd=lwd, col=linecol)
    for(i in 2:params$majordivs) {
      x1 <- x-adj[1]*wd+(i-1)*params$majordivplotunit
      y1 <- y-adj[2]*ht
      y2 <- y1+ht*tick.cex
      graphics::lines(c(x1, x1), c(y1, y2), col=linecol, lwd=lwd)
    }
  } else {
    stop("Invalid style specified to drawscalebar: ", style)
  }

}

#' Auto Plot Scalebar
#'
#' Automatically determines the geographical scale of the plot and
#' draws a labelled scalebar.
#'
#' @param plotunit The unit which the current plot is plotted in, one of \code{cm},
#' \code{m}, \code{km}, \code{in}, \code{ft}, \code{mi}. or \code{latlon}. This
#' parameter is optional if \code{plotepsg} is passed.
#' @param plotepsg The projection of the current plot. If extents are valid lat/lons,
#' the projection is assumed to be lat/lon (EPSG:4326), or Spherical Mercator otherwise
#' (EPSG:3857). This is done to work seamlessly with OpenStreetMap packages.
#' @param widthhint The fraction of the plottable width which the scale bar should
#' (mostly) occupy.
#' @param unitcategory One of "metric" or "imperial"
#' @param htin Height (in inches) of the desired scale bar
#' @param padin A vector of length 2 determining the distance in inches between the scalebar
#' and the edge of the plottable area.
#' @param style One of "bar" or "ticks".
#' @param bar.cols If \code{style=="bar"}, the colors to be repeated to make the bar.
#' @param lwd The line width to use when drawing the scalebar
#' @param linecol The line color to use when drawing the scalebar
#' @param tick.cex If \code{style=="ticks"}, the height of interior ticks.
#' @param labelpadin The distance between the end of the scalebar and the label (inches)
#' @param label.cex The font size of the label
#' @param label.col The color of the label
#' @param pos Where to align the scalebar. One of "bottomleft", "bottomright", "topleft",
#' or "topright".
#'
#' @export
#'
#' @examples
#' plot(1:5, 1:5, asp=1)
#' addscalebar(plotunit="m")
#' \donttest{
#' library(maptools)
#' data(wrld_simpl)
#' plot(wrld_simpl, xlim=c(-66.86, -59.75), ylim=c(43, 47.3)) #Nova Scotia
#' addscalebar()
#'
#' #also works in non-lat/lon coordinate systems
#' addscalebar(plotepsg=3395) #specify plot is in mercator projection
#' addscalebar(plotepsg=26920) #specify plot is in UTM Zone 20N
#'
#' }
#'
addscalebar <- function(plotunit=NULL, plotepsg=NULL, widthhint=0.25, unitcategory="metric",
                     htin=0.1, padin=c(0.15, 0.15), style="bar", bar.cols=c("black", "white"),
                     lwd=1, linecol="black", tick.cex=0.7, labelpadin=0.08, label.cex=0.8,
                     label.col="black", pos="bottomleft") {

  params <- scalebarparams(plotunit=plotunit, plotepsg=plotepsg, widthhint = widthhint,
                           unitcategory=unitcategory)
  extents <- params$extents

  bottomin <- graphics::grconvertY(extents[3], from="user", to="inches")
  leftin <- graphics::grconvertX(extents[1], from="user", to="inches")
  topin <- graphics::grconvertY(extents[4], from="user", to="inches")
  rightin <- graphics::grconvertX(extents[2], from="user", to="inches")

  ht <- graphics::grconvertY(bottomin+htin, from="inches", to="user") - extents[3]
  paduser <- graphics::grconvertX(leftin+labelpadin, from="inches", to="user") - extents[1]

  if(pos=="bottomleft") {
    x <- graphics::grconvertX(leftin+padin[1], from="inches", to="user")
    y <- graphics::grconvertY(bottomin+padin[2], from="inches", to="user")
    adj <- c(0,0)
    textadj <- c(0,0.5)
    textx <- x+params$widthplotunit+paduser
    texty <- y+0.5*ht
  } else if(pos=="topleft") {
    x <- graphics::grconvertX(leftin+padin[1], from="inches", to="user")
    y <- graphics::grconvertY(topin-padin[2], from="inches", to="user")
    adj <- c(0,1)
    textadj <- c(0, 0.5)
    textx <- x+params$widthplotunit+paduser
    texty <- y-0.5*ht
  } else if(pos=="topright") {
    x <- graphics::grconvertX(rightin-padin[1], from="inches", to="user")
    y <- graphics::grconvertY(topin-padin[2], from="inches", to="user")
    adj <- c(1,1)
    textadj <- c(1, 0.5)
    textx <- x-params$widthplotunit-paduser
    texty <- y-0.5*ht
  } else if(pos=="bottomright") {
    x <- graphics::grconvertX(rightin-padin[1], from="inches", to="user")
    y <- graphics::grconvertY(bottomin+padin[2], from="inches", to="user")
    adj <- c(1,0)
    textadj <- c(1, 0.5)
    textx <- x-params$widthplotunit-paduser
    texty <- y+0.5*ht
  }

  plotscalebar(x, y, ht, params, adj=adj, style=style, lwd=lwd, linecol=linecol,
               bar.cols=bar.cols, tick.cex=tick.cex)
  graphics::text(textx, texty, params$labeltext, adj=textadj, cex=label.cex, col=label.col)
}

