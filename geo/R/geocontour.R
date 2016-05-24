#' Plots contour lines.
#' 
#' The function plots contour lines. It is based on the Splus function
#' "contour" with some changes and extensions.  The main change is that
#' coordinates are in lat, lon instead of x, y. There are two additions : 1.
#' Possibility to have many colors.  2.  Possibility to only plot lines within
#' certain borders.  3.  Possibility to have labels in a box on the plot.
#' Note: to plot contour plots with fill between lines you use geocontour.fill
#' 
#' 
#' @param grd List with components \code{lat} and \code{lon} (or
#' \code{x},\code{y}) defining the grid values.  Can be made for example by lat
#' <- list(lat = seq(60,65,length = 100), lon = seq(-30,- 10,length = 200)).
#' Also lat can be the output of the program grid and then the components
#' \code{lat$grpt$lat}, \code{lat$grpt$lon} and \code{lat$reg} are used. See
#' geocontour.fill.
#' @param z Matrix or vector of values.  Length nlat*nlon except when lat is a
#' list with components \code{lat$grpt} and \code{lat$xgr}.  Then the length of
#' z is the same as the length of \code{lat$xgr$lat} and \code{lat$xgr$lon}.
#' @param nlevels Number of contourlines.  Used if the program has to determine
#' the contourlines.  Default value is 10. nlevels = 10 does not always mean 10
#' due to characteristic of the pretty command
#' @param levels Ratio of textsize to current size.  1 means same size, 0.5
#' half size , 0 no text etc.  Default value 1.
#' @param labcex Character expansion of letters indicating the z value on the
#' contour lines. Maybe used instead of a legend.
#' @param triangles If TRUE the program makes 4 triangles from each element in
#' the grid.  In fact it makes the number of grid points approximately twice as
#' many as before which means smoother contourlines.  Default value TRUE.
#' Triangles = TRUE is nessecary if the area was plotted with geocontour.fill.
#' It fits with geocontour.fill.
#' @param reg List with two components, \code{reg$lat} and \code{reg$lon}.
#' Region of interest. Contourlines are only plotted inside the region.  Holes
#' in the region of interest begin with NA.  If the region consists only of
#' holes then \code{reg$lat} and \code{reg$lon} begin with NA.
#' @param fill If fill is one the matrix is filled with zeros.  If two it is
#' filled with mean(z).  Default is fill = 1.  ( not an important parameter)
#' @param colors If colors is TRUE the contour lines are plotted in many
#' colors.  Else only in one.  The lines can also be in one color using many
#' linetypes.
#' @param col Number of color (colors) used to plot contourlines.  Default
#' value is 1 if colors is FALSE.  If colors is TRUE the default value is found
#' from the number of contour lines.
#' @param only.positive Logical value.  If FALSE then negative values are
#' allowed else negative values are set to zero.  Default value is FALSE.
#' @param maxcol Maximum color number in the Splus graphics, default value 155.
#' @param cex Character expansion of digits, default value is 0.7.
#' @param save If true the contour lines are saved in a list so they can later
#' be plotted with geolines.  Default value is false.
#' @param plotit If TRUE the lines are plotted on the graphics device else not.
#' Default value is true.
#' @param label.location List with components \code{lat} and \code{lon} (or
#' \code{x},\code{y}) Gives the lower left and upper right corner of the box
#' where the labels are put.  Default value is 0 that means no labels are put
#' on the drawing (except when \code{geopar$cont = TRUE}).  l1 is best given by
#' geolocator or directly by specifying label.location = "locator".
#' @param lwd Line with.  Default value is the value set when the program was
#' called.
#' @param lty ine type.  If lty is a vector of the same length of levels it
#' specifies the linetype for each contour line.  Default value is the same as
#' when the program was called.
#' @param labels.only If true only the labels are drawn.  Default is false.
#' @param digits Number of digits in labels.  Default value is one.
#' @param paint if true borders of regions will be painted.
#' @param set Set something.
#' @param col.names the names of the vectors containing the data in grd.
#' Default is col.names = c("lon","lat").
#' @param csi Size of character.  This parameter can not be set in R but for
#' compatibility with old Splus scripts the parameter cex is readjusted by cex
#' = cex*csi/0.12.  Use of this parameter is not recommended.  Default value is
#' NULL i.e not used.
#' @param drawlabels Draw labels on the contour lines? Default FALSE.
#' @return No values returned.
#' @section Side Effects: No side effects.
#' @seealso \code{\link{contour}}, \code{\link{geocontour.fill}},
#' \code{\link{geolocator}}, \code{\link{geopolygon}}, \code{\link{geotext}},
#' \code{\link{geosymbols}}, \code{\link{geogrid}}, \code{\link{geopar}},
#' \code{\link{geolines}}.
#' @examples
#' 
#'      ###################################################
#'      # Example l                                       #
#'      ###################################################  
#' \dontrun{
#'      geoplot(deg, cont = TRUE)                        # Plot initialized.
#'      geocontour(grd$lat,grd$lon,z,nlevels = 10,
#'                 neg = FALSE,reg = reg,colors = TRUE)          # Contour plot.
#'      geoplot(deg,pch = " ",cont = TRUE,new = TRUE)           # Plot over contourplot.
#' }
#'      ###################################################
#'      # Example 2 Sea Tempeture.                        #
#'      ###################################################  
#'       
#'      # The following data names used are in Icelandic, stodvar means
#'      # stations and botnhiti means temperature.
#' 
#'      geoplot()
#'      gbplot(500)
#'      grd <- list(lat = seq(63,67,length = 30),
#'                  lon = seq(-28,-10,length = 50))
#'      labloc <- list(lat = c(63.95,65.4),lon = c(-19.8,-17.3))
#' 
#'      grd1 <- geoexpand(grd)                       # Make grid.
#'      grd2 <- geoinside(grd1,gbdypif.500)  
#'      grd2 <- geoinside(grd2,island,robust = FALSE,option = 2) 
#'      # Use only the points where depth < 500 and outside Iceland.
#' \dontrun{
#'      #xx <- loess(botnhiti~ lat*lon,degree = 2,spaALSEn = 0.25,
#'      #            data = stodvar, na.action = na.omit)
#'      # Use loess for interpolating.
#' 
#'      #grd2$temp <- predict(xx,grd2)                
#'      #geocontour(grd2,z = grd2$temp,levels = c(0,1,2,3,4,5,6,7),
#'      #           label.location = labloc)
#' 
#'      ######################################################
#'      # Example 3 example of gam() and indexes.            #
#'      ###################################################### 
#' 
#'      stations<-data.frame(lat = stodvar$lat,lon = stodvar$lon,
#'                           temp = stodvar$botnhiti)
#'      # Making a partial data.frame from a big one called stodvar,
#'      # which means stations in Icelandic.
#' 
#'      stations[1:5,]             # Show first 5 lines all columns
#'                                 # in stations.
#'      dim(stations)              # Length of (lines,colums).
#'      dim(stations[!is.na(stations$temp),])      # Without NAs.
#'      my.data <- stations[!is.na(stations$temp),]
#'      my.data <- my.data[!is.na(my.data$lat),]
#'      my.data <- my.data[!is.na(my.data$lon),]
#'      # my.data is now same as stations but witout NAs in lat,
#'      # lon and temp.
#'        
#'      pred.grid <- list(lat = seq(63.25,67.25,length = round((67.25-63.25+1)*8)),
#'                        lon = seq(-27,-11.5,length = round((27-11.5+1)*4)))
#'      pred.grid <- geoexpand(pred.grid)
#'      # Making a grid to fit our area of interest.
#'      pred.grid <- geoinside(pred.grid,gbdypif.500)
#'      # Points within 500m.
#'      pred.grid <- geoinside(pred.grid,island,robust = FALSE,option = 2)
#'      # Points outside  of Iceland.
#'      
#'      geoplot(grid = FALSE)
#'      my.data <- geoinside(my.data,island,robust = FALSE,option = 2)
#'      geopoints(my.data,pch = ".")
#' 
#'      fit <- gam(temp~lo(lat,lon,span = 0.1),data = my.data)
#'      # see help(gam)
#'      # can also do:
#'      # fit <- loess(temp~lon*lat,data = my.data,span = 0.1)
#'      # fit <- gam(temp~ns(lon,df = 7)*ns(lat,df = 5),data = my.data) 
#' 
#'      pred.grid$pred.temp <- predict(fit,newdata = pred.grid)
#'      geocontour(pred.grid,z = pred.grid$pred.temp,levels = 0:7,
#'                 label.location = labloc)
#' }
#' @export geocontour
geocontour <-
function (grd, z, nlevels = 10, levels = NULL, labcex = 1, triangles = TRUE, 
    reg = 0, fill = 1, colors = TRUE, col = 1, only.positive = FALSE, 
    maxcol = 155, cex = 0.7, save = FALSE, plotit = TRUE, label.location = 0, 
    lwd = 1, lty = 1, labels.only = FALSE, digits = 1, paint = FALSE, 
    set = NA, col.names = c("lon", "lat"), csi = NULL, drawlabels = FALSE)
{
    geopar <- getOption("geopar")
    if (!is.null(csi)) 
        cex <- cex * csi/0.12
    if (!is.null(attributes(grd)$grid)) {
        z <- grd
        grd <- attributes(grd)$grid
    }
    limits <- NULL
    maxn <- 10000
    grd <- Set.grd.and.z(grd, z, NULL, set, col.names)
    z <- grd$z
    z <- z + rnorm(length(z)) * 1e-09
    grd <- grd$grd
    grd <- extract(grd, z, maxn, limits, col.names = col.names)
    z <- grd$z
    grd <- grd$grd1
    ind <- c(1:length(z))
    ind <- ind[is.na(z)]
    if (length(ind) > 0) {
        if (fill == 0) 
            z[ind] <- -99999
        if (fill == 1) 
            z[ind] <- 0
        if (fill == 2) 
            z[ind] <- mean(z)
    }
    lon <- grd[[col.names[1]]]
    lat <- grd[[col.names[2]]]
    if (only.positive) {
        ind <- c(1:length(z))
        ind <- ind[z < mean(z[z > 0])/1000 & z != -99999]
        z[ind] <- mean(z[z > 0])/1000
    }
    cond1 <- col.names[1] == "lon" && col.names[2] == "lat"
    cond2 <- col.names[1] == "x" && col.names[2] == "y" && geopar$projection == 
        "none"
    if (cond1 || cond2) {
        oldpar <- selectedpar()
        on.exit(par(oldpar))
        par(geopar$gpar)
        if (geopar$cont) 
            par(plt = geopar$contlines)
    }
    if (cex != 0) 
        par(cex = cex)
    nx <- length(lon)
    ny <- length(lat)
    lon1 <- matrix(lon, nx, ny)
    lat1 <- t(matrix(lat, ny, nx))
    if (!labels.only) {
      if (geopar$projection == "Mercator" && col.names[1] == 
          "lon") {
        z <- matrix(z, nrow = length(lon), ncol = length(lat))
        lon2 <- c(matrix(lon[1], length(lat), 1))
        lat2 <- c(matrix(lat[1], length(lon), 1))
        xlat <- Proj(lat, lon2, geopar$scale, geopar$b0, 
                     geopar$b1, geopar$l1, geopar$projection)
        xlon <- Proj(lat2, lon, geopar$scale, geopar$b0, 
                     geopar$b1, geopar$l1, geopar$projection)
      }
      else {
        z <- matrix(z, nrow = length(lon), ncol = length(lat))
        xlon <- list(x = lon)
        xlat <- list(y = lat)
      }
    }
    if (colors) {
      if (is.null(levels)) {
        if (nlevels == 0) 
          nlevels <- 10
        levels <- pretty(z, nlevels)
      }
      nlevels <- length(levels)
      if (length(lty) == length(levels) && length(levels) > 
          1) 
        linetypes <- TRUE
      else {linetypes <- FALSE;lty <- rep(lty,length(levels))}
      if (length(lwd) == length(levels) && length(levels) > 
          1) 
        linew <- TRUE
      else {linew <- FALSE;lwd <- rep(lwd,length(levels))}
      if (length(col) == 1) {
        if (length(lty) == length(levels)) 
          color <- rep(1, nlevels)
        else {
          mincol <- 2
          color <- c(1:nlevels)
          color <- round(2 + ((color - 1) * maxcol)/(nlevels))
        }
      }
      else color <- col
      if (!labels.only) {
        if (length(ind) > 1) 
          z[ind] <- NA
        if (geopar$projection == "Lambert") {
          lev <- contourLines(lon + 400, lat, z, levels = levels)
          for (i in 1:length(lev)) {
            j <- match(lev[[i]]$level,levels)
            if (linew) 
              lw <- lwd[j]
            else lw <- 0
            if (linetypes) 
              lt <- lty[j]
            else lt <- 0
            geolines(lev[[i]]$y, lev[[i]]$x - 400, col = color[j], 
                     lty = lt, lwd = lw)
          }
        }
        else {
          for (i in 1:nlevels) {
            lev <- contour(xlon$x, xlat$y, z, axes = FALSE, drawlabels=drawlabels,
                           levels = c(levels[i], levels[i]), add = TRUE, 
                           triangles = triangles, labcex = labcex, xlim = geopar$limx, 
                           ylim = geopar$limy, col = color[i], xlab = " ", 
                           ylab = " ", save = FALSE, plotit = TRUE,lwd=lwd[i],lty=lty[i])
          }
        }
      }
    }
    else {
      if (length(ind) > 1) 
        z[ind] <- NA
      if (geopar$projection == "Lambert") {
        if (length(levels) == 1) 
          lev <- contourLines(lon + 400, lat, z, nlevels = nlevels)
        else {
          lev <- contourLines(lon + 400, lat, z, levels = levels)
          for (i in 1:length(lev)) {
            geolines(lev[[i]]$y, lev[[i]]$x - 400, col = col)
          }
          
        }
      }
      else {
        if (length(levels) == 1) 
          lev <- contour(xlon$x, xlat$y, z, nlevels = nlevels,
                         triangles = triangles, labcex = labcex, add = TRUE, 
                         xlim = geopar$limx, ylim = geopar$limy, axes = FALSE, 
                         col = col, xlab = " ", ylab = " ", save = save, 
                         plotit = plotit ,drawlabels=drawlabels)
        else lev <- contour(xlon$x, xlat$y, z, axes = FALSE, 
                            levels = levels, triangles = triangles, add = TRUE, 
                            labcex = labcex, xlim = geopar$limx, ylim = geopar$limy, 
                            col = col, xlab = " ", ylab = " ", save = save, 
                            plotit = plotit, drawlabels=drawlabels)
      }
    }
    if (save && geopar$projection == "Mercator") {
      lev <- contourLines(xlon$x, xlat$y, z, levels = levels)
      tmpdata <- data.frame(lat=0,lon=0,level=-99)
      res <- NULL
      for (i in 1:length(lev)) {
        tmp <- invProj(lev[[i]])
        tmp <- data.frame(lat=tmp$lat,lon=tmp$lon)
        tmp$level <- rep(lev[[i]]$level, nrow(tmp))
        res <- rbind(res, tmp)
        res <- rbind(res,tmpdata)
      }
      i <- res$lat == 0
      if(any(i)) res$lat[i] <- res$lon[i] <- NA
      lev <- res
    }
    if (save && geopar$projection == "Lambert") {
      tmpdata <- data.frame(lat=0,lon=0,level=-99)
      res <- NULL
      for (i in 1:length(lev)) {
        tmp <- lev[[i]]
        tmp <- data.frame(lat=tmp$y,lon=tmp$x-400)
        tmp$level <- rep(lev[[i]]$level, nrow(tmp))
        res <- rbind(res, tmp)
        res <- rbind(res,tmpdata)
      }
      i <- res$lat == 0
      if(any(i)) res$lat[i] <- res$lon[i] <- NA
      lev <- res
    }

    if (geopar$projection == "Lambert") 
      par(geopar$gpar)
    if (length(label.location) == 1) 
      if (label.location == "locator") 
        label.location <- geolocator(n = 2)
    if (length(label.location) > 1) {
      label.location <- Proj(label.location, scale = geopar$scale, 
                             b0 = geopar$b0, b1 = geopar$b1, l1 = geopar$l1, projection = geopar$projection)
      if (geopar$projection == "none") 
        paint.window.x(label.location, border = TRUE)
      else paint.window(label.location, border = TRUE)
        labels_line(levels, digits, color, lty, lwd, xlim = label.location$x, 
                    ylim = label.location$y, linew)
    }
    if (geopar$cont && colors) {
      par(plt = geopar$contlab)
      par(new = TRUE)
      plot(c(0, 1, 1, 0, 0), c(0, 0, 1, 1, 0), type = "l", 
           axes = FALSE, xlab = " ", ylab = " ")
      labels_line(levels, digits, color, lty, linew)
    }
    if (length(reg) > 1 && paint) {
      nx <- length(lon)
      ny <- length(lat)
      lon <- matrix(lon, nx, ny)
      lat <- t(matrix(lat, ny, nx))
      shadeborder(reg, lat, lon)
    }
    if (cond1 || cond2) {
      par(oldpar)
    }
    if (save) {
      return(invisible(lev))
    }
    else return(invisible())
  }

