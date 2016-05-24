ColorMap <- structure(function#Plot Levels of a Variable in a Colour-Coded Map
### Plot Levels of a Variable in a Colour-Coded Map
(
  values, ##<< variable to plot
  map=NULL, ##<< map object
  polys=NULL, ##<< an object of class SpatialPolygons (See \link[sp]{SpatialPolygons-class}
  log = FALSE, ##<< boolean of whether to plot values on log scale
  nclr = 7, ##<< number of colour-levels to use
  include.legend = list(TRUE), ##<< boolean of whether to include legend
  round = 3, ##<< number of digits to round to in legend
  brks = NULL, ##<< if desired, pre-specified breaks for legend
  legend = NULL, ##<< if desired, a pre-specified legend
  location = "topright", ##<< location of legend
  rev = FALSE, ##<< boolean of whether to reverse colour scheme (darker colours for smaller values)
  alpha = 0.5, ##<< alpha value of colors
  GRAY = FALSE, ##<< boolean: if TRUE, use gray scale instead
  palette = c("YlOrRd", "RdYlGn","Spectral")[1],##<< palette to choose from RColorBrewer
  textInPolys = NULL, ##<< text to be displayed inside polygons. This can be a column names for values
  ... ##<< extra args to pass to \code{PlotPolysOnStaticMap}
  ){
    if (length(palette) == nclr) {
      plotclr <- palette
    } else if (GRAY | !requireNamespace("RColorBrewer", quietly = TRUE) ) {
      plotclr <- grey(1 - seq(0, 1, by = 1/(nclr - 1)))
    } else {
      plotclr <- RColorBrewer::brewer.pal(nclr,palette)
      #display.brewer.all()
    }
   
    plotclr = AddAlpha(plotclr,alpha)
    
    nclr <- nclr + 1
    if (is.null(brks)) {
        if (log) {
            brks <- exp(seq(from = min(log(values)), to = max(log(values)), 
                length.out = nclr))
        }
        else {
            brks <- seq(from = min(values), to = max(values), 
                length.out = nclr)
        }
    }
    nclr <- nclr - 1
    print(brks)
    if (rev) {
        plotclr <- rev(plotclr)
    }
    colornum <- findInterval(values, brks, all.inside = T)
    colcode <- plotclr[colornum]
    if (is.null(legend)) {
      leglabs =
      function (vec, under = "under", over = "over", between = "-") 
      {
        x <- vec
        lx <- length(x)
        if (lx < 3) 
          stop("vector too short")
        res <- character(lx - 1)
        res[1] <- paste(under, x[2])
        for (i in 2:(lx - 2)) res[i] <- paste(x[i], between, x[i + 
          1])
        res[lx - 1] <- paste(over, x[lx - 1])
        res
      }
        legend <- leglabs(signif(brks, digits = round))
        legend[1] <- paste(signif(min(values), digits = round), 
            "-", substr(legend[1], 7, nchar(legend[1])))
        legend[nclr] <- paste(substr(legend[nclr], 6, nchar(legend[nclr])), 
            "-", signif(max(values), digits = round))
    }
    if (is.null(polys)) {#no plotting, just return the colors
      return(list(colcode=colcode, legend = legend, fill = plotclr))
    } else {
      if (is.null(map)) {
        PBSmapping::plotPolys(polys, col = colcode )#hoping, this is well defined by the class somehow !
      } else {
        if (inherits(polys, 'Spatial')) polys = SpatialToPBS(polys)$xy
        #if (!is.null(textInPolys)) textInPolys = as.character(values[,textInPolys])
        PlotPolysOnStaticMap(map, polys, col = colcode, textInPolys=textInPolys,...)
      }
    
    if (include.legend[[1]]) 
        legend(location, legend = legend, fill = plotclr, bty = "n")
    }
}, ex = function(){
  if (interactive()){
  data("NYleukemia", envir = environment())
  population <- NYleukemia$data$population
  cases <- NYleukemia$data$cases
  mapNY <- GetMap(center=c(lat=42.67456,lon=-76.00365), destfile = "NYstate.png", 
                  maptype = "mobile", zoom=9)
  ColorMap(100*cases/population, mapNY, NYleukemia$spatial.polygon, add = FALSE,
           alpha = 0.35, log = TRUE, location = "topleft")
}
  #ColorMap(100*cases/population, map=NULL, NYleukemia$spatial.polygon)
  
})
