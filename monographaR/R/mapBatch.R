mapBatch <-
function (data, zoom = T, margin = 0.1, axes = T, shape=NULL, 
          export = "pdf", raster = NULL, points.col = "black", 
          points.border = "gray50", points.cex = 1, shape.col = "white", 
          shape.border = "black", raster.col = rev(gray.colors(65, start = 0, end = 1)),
          raster.legend=F, hillshade = F, width = 8, height = 8, image.resolution = 100, 
          figure.number = T, title = T, box = T, add.minimap = F, minimap.shape=NULL,
          minimap.shape.col = "white", minimap.shape.border = "gray50", 
          minimap.pos = "topleft", minimap.add.points = T, minimap.points.col = "black",
          minimap.points.border = "gray50", minimap.points.cex = 1, minimap.extent = NULL,
          maxpixels=100000, ...) 
{
  if (class(data) != "data.frame") {
    stop("data must be a data.frame")
  }
  if (ncol(data) != 3) {
    stop("data must have 3 columns, see help(\"mapBatch\")")
  }
  wrld_simpl = NULL	
  if (is.null(shape)) {
    data(wrld_simpl, envir = environment())
    wrld_simpl -> shape
  }
  if (add.minimap) { 
    try(dev.off(), silent=T)
    zoom = T 
  }
  if (is.null(raster) == F && hillshade == T) {
    proj4string(raster) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
    slope <- terrain(raster, opt='slope')
    aspect <- terrain(raster, opt='aspect')
    hill <- hillShade(raster, aspect, 40, 270)
    rev(gray.colors(100, start = 0, end = 1, alpha = 0.3)) -> hill.col
  }
  message("Assuming the columns are ordered as: species, longitude and latitude")
  colnames(data) <- c("sp", "x", "y")
  geo <- data
  coordinates(geo) <- ~x + y
  min.margin = 5
  max.margin = 10
  mmm <- (min.margin++(margin/mean(c(extent(geo)[2]-extent(geo)[1], extent(geo)[4]-extent(geo)[3]))))
  if (mmm > max.margin) {
    mmm <- max.margin
  }
  ext.all <- extent(geo)++mmm
  spp <- as.character(unique(data[, 1]))
  spp <- sort(spp)
  if (export == "pdf") {
    pdf("mapBatch.pdf", width=width, height=height)
  }
  for (i in 1:length(spp)) {
    sp <- spp[i]
    if (export == "tiff") {
      if (figure.number) {
        paste("Figure ", i, " - ", sp, ".tif", sep = "") -> lab0
      } else {
        paste(sp, ".tif", sep="") -> lab0
      }
      tiff(lab0, width = width, height = height, units = "in", 
           res = image.resolution, compression = "lzw+p")
    }
    if (export == "jpeg") {
      if (figure.number) {
        paste("Figure ", i, " - ", sp, ".jpg", sep = "") -> lab0
      } else {
        paste(sp, ".jpg", sep="") -> lab0
      }
      jpeg(lab0,width = width, height = height, units = "in", 
          quality = 90, res = image.resolution)
    }
    spRows <- which(data$sp == sp)
    spData <- data[spRows, ]
    xy <- spData
    coordinates(xy) <- ~x + y
    if (zoom == T) {
      ext <- extent(xy)++mmm
      xlim <- c(ext[1], ext[2])
      ylim <- c(ext[3], ext[4])
    }
    else {
      ext <- ext.all
      xlim <- c(ext[1], ext[2])
      ylim <- c(ext[3], ext[4])
    }
    if (class(shape) == "list") {
      plot(shape[[1]], xlim = xlim, ylim = ylim, axes = axes, 
           col = shape.col, border = shape.border, add = F, 
           asp = 1, ...)
      if (is.null(raster) == F) {
        if (hillshade == T) {
          plot(hill, col=hill.col, legend=F, axes=F, box=F, add=T, maxpixels=maxpixels)
          plot(raster, col = raster.col, add = T, legend=raster.legend, maxpixels=maxpixels)
        } else {
          plot(raster, col = raster.col, add = T, legend=raster.legend, maxpixels=maxpixels)
        }
        plot(shape[[1]], xlim = xlim, ylim = ylim, axes = axes, 
             col = shape.col, border = shape.border, add = T, 
             asp = 1, ...)
      }
      for (k in 2:length(shape)) {
        plot(shape[[k]], xlim = xlim, ylim = ylim, 
             axes = axes, col = shape.col, border = shape.border, 
             add = T, asp = 1)
      }
    } else {
      plot(shape, xlim = xlim, ylim = ylim, axes = axes, 
           col = shape.col, border = shape.border, add = F, 
           asp = 1)
      if (is.null(raster) == F) {
        if (hillshade == T) {
          plot(hill, col=hill.col, legend=F, axes=F, box=F, add=T, maxpixels=maxpixels)
          plot(raster, col = raster.col, add = T, legend=raster.legend, maxpixels=maxpixels)
        } else {
          plot(raster, col = raster.col, add = T, legend=raster.legend, maxpixels=maxpixels)
        }
        plot(shape, xlim = xlim, ylim = ylim, axes = axes, 
             col = shape.col, border = shape.border, add = T, 
             asp = 1)
      }
    }
    plot(xy, pch = 21, col = points.border, bg = points.col, cex = points.cex, add = T)
    if (box) { box() }
    if (title) { title(sp) }
    if (add.minimap) {
      if (is.null(minimap.shape)) { 
        if (class(shape) == "list") {
          shape[[length(shape)]] -> minimap.shape
        } else {
          shape -> minimap.shape
        }
      }
      par()$usr -> lims
      mw = min(width,height)*0.2
      if (is.null(minimap.extent)) {
        extent(geo)++1 -> ext.m
      } else { minimap.extent -> ext.m }
      xlim.m <- c(min(ext.m[1], lims[1]), max(ext.m[2], lims[2]))
      ylim.m <- c(min(ext.m[3], lims[3]), max(ext.m[4], lims[4]))
      png("temp.png", bg='transparent', width=mw, height=mw, units="in", res=400)
      par(mar=c(0,0,0,0))
      plot(minimap.shape, xlim=xlim.m, ylim=ylim.m, col=minimap.shape.col, border=minimap.shape.border)
      if (minimap.add.points) {
        plot(xy, pch = 21, col = minimap.points.border, bg = minimap.points.col, cex = minimap.points.cex, add = T)
      }
      rect(lims[1], lims[3], lims[2], lims[4], lwd=2, lty="dotted")
      dev.off()
      jpg0 = readPNG("temp.png", native=T)
      min(c(lims[2]-lims[1], lims[4]-lims[3])) -> min.plot
      if (minimap.pos == "topleft") {
        lims[1]++(min.plot*0.02) -> x1
        lims[1]++(min.plot*0.27) -> x2
        lims[4]-(min.plot*0.02) -> y2
        lims[4]-(min.plot*0.27) -> y1
      }
      if (minimap.pos == "bottomright") {
        lims[2]-(min.plot*0.27) -> x1
        lims[2]-(min.plot*0.02) -> x2
        lims[3]++(min.plot*0.27) -> y2
        lims[3]++(min.plot*0.02) -> y1
      }
      rasterImage(jpg0, x1, y1, x2, y2, angle=0, interpolate=T)
      unlink("temp.png")
    }
    if (export == "tiff" || export == "jpeg") { dev.off() }
  }
  if (export == "pdf") { dev.off() }
  cat("Maps were saved in:")
  cat("\n", getwd())
}
