PlotStereoMap <- function(var, lon, lat, latlims = c(60,90), toptitle = '',  
                          sizetit = 1, units = '', brks = NULL, cols = NULL,
                          filled.continents = FALSE, dots = NULL, intlat=10, 
                          drawleg = TRUE, subsampleg = 1, colNA='white') {
  #
  #  Input arguments 
  # ~~~~~~~~~~~~~~~~~
  #
  dims <- dim(var)
  if (length(dims) > 2) {
    stop("Only 2 dimensions expected for var : (lon,lat) ")
  }
  if (dims[1] != length(lon) | dims[2] != length(lat)) {
    if (dims[1] == length(lat) & dims[2] == length(lon)) {
      var <- t(var)
      dims <- dim(var)
    } else {
      stop("Inconsistent var dimensions / longitudes +  latitudes")
    }
  }
  dlon <- median(lon[2:dims[1]] - lon[1:(dims[1] - 1)])/2
  dlat <- median(lat[2:dims[2]] - lat[1:(dims[2] - 1)])/2
  colorbar <- colorRampPalette(c("dodgerblue4", "dodgerblue1", 
              "forestgreen", "yellowgreen", "white", "white", 
              "yellow", "orange", "red", "saddlebrown"))
  if (is.null(brks) == T) {
    ll <- signif(min(var, na.rm = T), 4)
    ul <- signif(max(var, na.rm = T), 4)
    if (is.null(cols) == T) {
      cols <- colorbar(10)
    }
    nlev <- length(cols)
    brks <- signif(seq(ll, ul, (ul - ll) / nlev), 4)
  } else {
    if (is.null(cols) == T) {
      nlev <- length(brks) - 1
      cols <- colorbar(nlev)
    } else {
      if (length(cols) != (length(brks) - 1)) {
        stop("Inconsistent colour levels / list of colours")
      }
    }
  }
  #
  #  Defining the layout
  # ~~~~~~~~~~~~~~~~~~~~~
  # 
  xmargin <- 0
  ymargin <- 0
  topmargin <- 0
  ymargin2 <- 0 
  if (toptitle != '') {
    topmargin <- 1.5 
  }
  #
  #  Plotting the map
  # ~~~~~~~~~~~~~~~~~~
  #
  if (drawleg) {
    layout(matrix(1:2, ncol = 2, nrow = 1), widths = c(9, 1))
  }
  par(mar = c(xmargin, ymargin, topmargin, ymargin2), mgp = c(1, 1, 0) ,cex = 1.4)
  map("world", projection = "stereographic", orientation = c(sign(latlims[1])*90, 0, 0),
       xlim = c(-180,180), ylim = latlims, interior = F, wrap = T, lwd = 1)
  for (jx in 1:dim(var)[1]) { 
    for (jy in 1:dim(var)[2]) {
      if ( lat[jy] >= latlims[1] & latlims[2] >= lat[jy] ) {
        coord=mapproject(c(lon[jx]-dlon,lon[jx]+dlon,lon[jx]+dlon,lon[jx]-dlon),
                         c(lat[jy]-dlat,lat[jy]-dlat,lat[jy]+dlat,lat[jy]+dlat))
        ind=which(brks[-1]>var[jx,jy] & var[jx,jy]>brks[-length(brks)])
        if (is.na(var[jx,jy])) {
          polygon(coord, col=colNA, border=NA)
        }else{
          polygon(coord, col=cols[ind], border=NA)
        }
        if (is.null(dots) == F) {
          if (dots[jx, jy] == T) {
            coord=mapproject(lon[jx],lat[jy])
            points(coord, col = 'black', cex = 0.02, pch = 4)
          }
        }
      }
    }
  }
  if (filled.continents == T) {
    map("world", projection = "stereographic", orientation = c(sign(latlims[1])*90, 0, 0), 
        xlim = c(-180,180), ylim = latlims, interior = F, wrap = T, lwd = 0.7, 
        col = gray(0.5), fill = T, add = T, border = gray(0.5))
  }else{
    map("world", projection = "stereographic", orientation = c(sign(latlims[1])*90, 0, 0),
       xlim = c(-180,180), ylim = latlims, interior = F, wrap = T, lwd = 1, add = T)
  }
  map.grid(lim = c(-180, 180, latlims), nx = 18, ny = (latlims[2]-latlims[1])/intlat,
           col = 'lightgrey', labels = FALSE) 
  title(toptitle, cex.main=sizetit)
  #
  #  Colorbar
  # ~~~~~~~~~~
  #
  if (drawleg) {
    par(mar = c(xmargin+0.3, 0, topmargin, 2.5), mgp = c(1, 1, 0),
        las = 1, cex = 1.2)
    image(1, c(1:length(cols)), t(c(1:length(cols))), axes = F, col = cols,
          xlab = '', ylab = '', main = units, cex.main = 1.1)
    box()
    axis(4, at = seq(0.5, length(brks) - 0.5, subsampleg),
         labels = brks[seq(1, length(brks), subsampleg)])
  }
}




