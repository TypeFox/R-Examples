## Labels drawing
## TODO: labels' rotations.
## first, in no boxes, it is easy
## if boxes, at least do  90 degrees rotations
## finally, more than one rotation possible.
adeg.panel.label <- function(x, y, labels, plabels, pos = NULL) {
  if(any(plabels$cex > 0)) {
    n <- length(x)
    plboxes <- plabels$boxes  
    draw <- plabels$cex > 0
		## using .textsize funtion in utils.R
    textS <- .textsize(labels, plabels)
    srt <- textS$srt
    if(plboxes$draw && (srt != 0) & srt != 90) 
      warning("Boxes only implemented for 0 or 90 degrees rotation", call. = FALSE)
    ldraw <- rep(draw, length.out = n)  ## draw long enough
    ldraw[which(is.na(labels[1:n]) | labels[1:n] == "")] <- FALSE ## if no labels or null string don't bother
    width <- rep(textS$w, length.out = n)[ldraw]
    height <- rep(textS$h, lenght.out = n)[ldraw]
    lab <- rep(labels, length.out = n)[ldraw] ## no NA, removed using ldraw
    bdraw <- rep(plboxes$draw, length.out = length(ldraw))
    
    ## no boxes if no labels
    bdraw <- (bdraw & ldraw)
    ## labels a dessiner
    optim <- plabels$optim[1] ## only one possibility
    newpos <- list(x = x[ldraw], y = y[ldraw])
    
    if(optim) {                
      ## calcul des nouvelles positions uniquement pour les labels qui seront dessines
      ## informations sur panel
      nativelim <- current.panel.limits(unit = "native")
      incheslim <- current.panel.limits(unit = "inches")
      ## calcul des nouvelles positions.
      if(any(is.na(width)) | any(is.na(height)) | any(is.na(newpos$y)) | any(is.na(newpos$x))) 
        stop("NA restants revoir adeg.panel.label")
      newpos <- .pointLabel(x = newpos$x, y = newpos$y, labels = lab, width = width / diff(nativelim$xlim), height = height / diff(nativelim$ylim), limits = nativelim, xyAspect = diff(incheslim$xlim) / diff(incheslim$ylim), trace = FALSE)
    }
    
    if(any(bdraw)) {
      ## dessins de chaque boite avec son label
      plboxes <- lapply(plboxes, FUN = function(x) {rep(x, length.out = length(ldraw))})
      srt <- rep(srt, length.out = length(ldraw))
      plabels <- lapply(plabels, FUN = function(x) {rep(x, length.out = n)[ldraw]})
      for(i in 1:length(newpos$x)) {
        if(bdraw[ldraw][i]) {
          ## labels sizes
          panel.rect(
            x = unit(newpos$x[i], "native"),
            y = unit(newpos$y[i], "native"),
            width = width[i],
            height = height[i],
            col = plboxes$col[ldraw][i],
            alpha = plboxes$alpha[ldraw][i],
            border = plboxes$border[ldraw][i],
            lty = plboxes$lty[ldraw][i],
            lwd = plboxes$lwd[ldraw][i]
          )
        }
        panel.text(labels = lab[i], x = unit(newpos$x[i], "native"), y = unit(newpos$y[i], "native"), col = plabels$col[i], cex = plabels$cex[i], alpha = plabels$alpha[i], srt = srt[i])
      }
    }
    else { ## only text
      if(any(!ldraw)) ## obliger de repeter pour dessiner si un label doit etre ignorer
        panel.text(labels = lab, x = unit(newpos$x, "native"), y = unit(newpos$y, "native"),
          				 col = rep(plabels$col, length.out = length(ldraw))[ldraw], cex = rep(plabels$cex, length.out = length(ldraw))[ldraw], 
          				 alpha = rep(plabels$alpha, length.out = length(ldraw))[ldraw], rep(srt, length.out = length(ldraw))[ldraw], pos = pos)
      else
        panel.text(labels = lab, x = unit(newpos$x, "native"), y = unit(newpos$y, "native"), 
          				 col = plabels$col, cex = plabels$cex, alpha = plabels$alpha, srt = srt, pos = pos)
    }
  }
}


adeg.panel.nb <- function(nbobject, coords, col.edge = "black", lwd = 1, lty = 1, pch = 20, cex = 1, col.node = "black", alpha = 1) {
  if(class(nbobject) != "nb")
    stop("nb object is not class nb") ## prevoir dans les fonctions user une selection de l element neighbourght si object de type listw
  if(length(nbobject) != nrow(coords))
    stop("error for nb object, not the same numbers of nodes and coordinates", call. = FALSE)
  edges <- cbind(rep(1:length(nbobject), lapply(nbobject, length)), unlist(nbobject))
  edges <- edges[edges[,2] != 0, ]
  
  ## ici faire rep des parametres pour pouvoir ensuite modifier couleur
  adeg.panel.edges(edges, coords, col.edge, lwd, lty, pch, cex, col.node, alpha)
}


## adeg.panel.edges....
## col, lwd, lty etc peuvent varier selon poids des connexions
adeg.panel.edges <- function(edges, coords, col.edge = "black", lwd = 1, lty = 1, pch = 20, cex = 1, col.node = "black", alpha = 1) {
  panel.points(x = coords[, 1], y = coords[, 2], col = col.node, pch = pch, alpha = alpha, cex = cex)
  panel.segments(x0 = coords[edges[, 1], 1], y0 = coords[edges[, 1], 2], x1 = coords[edges[, 2], 1], y1 = coords[edges[, 2], 2], col = col.edge, lwd = lwd, lty = lty)
}


################## Panel.spatial #############################
## spObject can be :
## SpatialGridDataFrame","SpatialLinesDataFrame","SpatialPixelsDataFrame","SpatialPointsDataFrame","SpatialPolygonsDataFrame"
## n : nombre intervales si data 
## TODO: spObject pourrait etre une liste
adeg.panel.Spatial <- function(SpObject, sp.layout = NULL, col = 1, border = 1, lwd = 1, lty = 1, alpha = 0.8, cex = 1, pch = 20, n = length(col), spIndex = 1, ...) {

  if(length(grep("DataFrame", class(SpObject))) > 0) { ## there is data in 'SpObject' (it is a SpatialPolygonsDataFrame).
    mapSp <- try(SpObject[names(SpObject)[spIndex]], silent = TRUE) ## only the first map (spIndex = 1)
    values <- try(mapSp@data[, 1], silent = TRUE)
    
    if(is.factor(values)) { ## qualitative values
      if(length(col) != nlevels(values)) {
        if(length(col) == 1)  ## all values have the same color
          col <- rep(col, length.out = nlevels(values))
        else 
        	col <- adegpar()$ppalette$quali(nlevels(values))
      	colvalue <- col[values]
      } else
        colvalue <- col
    
    } else {  ## quantitative values
      breaks <- pretty(values, length(col))
      if((length(breaks) - 1) != length(col)) {
        if(length(col) == 1)  ## 'col' is not modified by the user
        	col <- adegpar()$ppalette$quanti(length(breaks) - 1)
        else  ## 'col' is modified but there is not enough color values
          col <- colorRampPalette(col)(length(breaks) - 1)
      }
      colvalue <- col[cut(values, breaks, include.lowest = TRUE)]
    }
  } else {  ## there is no data in 'SpObject'
    mapSp <- SpObject
    colvalue <- col
  }

   if(!is.null(sp.layout))
     sppanel(sp.layout)
  
  if(inherits(SpObject, what = "SpatialPoints")) {
    ## insert ppoints.parameters for pch and cex
    sp.points(mapSp, col = colvalue, pch = pch, cex = cex, alpha = alpha)
  }
  
  if(inherits(SpObject, what = "SpatialPolygons"))
    sp.polygons(mapSp, col = border, fill = colvalue, alpha = alpha, lty = lty, lwd = lwd)
  
  ## For spatialLine problems ; no various colors
  if(inherits(SpObject, what = "SpatialLines"))
    sp.lines(mapSp, col = colvalue, alpha = alpha, lty = lty, lwd = lwd)
  if(inherits(SpObject, what = "SpatialGrid"))
    sp.grid(mapSp, at = breaks, col = col)
}


adeg.panel.values <- function(x, y, z, method, symbol, ppoints, breaks, centerpar = NULL, center = 0) {
  if((length(x) != length(y)) | (length(y) != length(z)))
    stop("error in panel.values, not equal length for x, y, and z")
  
  maxsize <- max(abs(breaks))  ## biggest value
  z <- z - center

  if(!missing(center) & !is.null(centerpar)) {
      xnull <- x[abs(z) < sqrt(.Machine$double.eps)]
      ynull <- y[abs(z) < sqrt(.Machine$double.eps)]
  }
  
  if(method == "size"){
      size <- .proportional_map(z, maxsize) * ppoints$cex[1]
      colfill <- ifelse(z < 0, ppoints$col[1], ppoints$col[2])
      colborder <- ifelse(z < 0, ppoints$col[2], ppoints$col[1])
      
  } else if(method == "color"){
      size <- ppoints$cex[1]
      breaks <- sort(breaks)
      colfill <- ppoints$col[as.numeric(cut(z, breaks, include.lowest = TRUE))]
      if(any(is.null(colfill)) | any(is.na(colfill)))
             stop("error in the definition of color symbol", call. = FALSE)
      colborder <- "black"
  }
  
  cstnormal <- 5 ## same value in createkey
  panel.points(x = x, y = y, cex = size * cstnormal, pch = .symbol2pch(symbol), fill = colfill, col = colborder, alpha = ppoints$alpha)
  if(!missing(center) && !is.null(centerpar))
    panel.points(x = xnull, y = ynull, pch = centerpar$pch, col = centerpar$col, cex = centerpar$cex)
  return(cstnormal) 
}


adeg.panel.hist <- function(histValues, horizontal = TRUE, densi, drawLines, params = list(), identifier = "histogramADEg") {
  ## from panel.histogram of the lattice package
  plot.polygon <- modifyList(list(plot.polygon = trellis.par.get("plot.polygon")), params, keep.null = TRUE)[[1L]] ## hist params
  add.line <- modifyList(list(add.line = trellis.par.get("add.line")), params, keep.null = TRUE)[[1L]] ## backgroundlines
  plot.line <- modifyList(list(plot.line = trellis.par.get("plot.line")), params, keep.null = TRUE)[[1L]] ## density line
  
  h <- histValues
  breaks <- h$breaks
  heiBar <- h$counts
  nb <- length(breaks)
  ## counts lines
  if(horizontal)
    do.call("panel.abline", c(list(h = drawLines), add.line))
  else
    do.call("panel.abline", c(list(v = drawLines), add.line))
  ## warning : density lines re-scale to check
  contdensi <- (h$counts[h$density != 0 & h$counts != 0] / h$density[h$density != 0 & h$counts != 0])[1]
  if(horizontal) {
    if(nb > 1) {
      panel.rect(x = h$mids, y = 0, height = heiBar, width = diff(breaks), 
                 col = plot.polygon$col, alpha = plot.polygon$alpha, border = plot.polygon$border, lty = plot.polygon$lty, 
                 lwd = plot.polygon$lwd, just = c("center", "bottom"), identifier = identifier)
    }
    do.call("panel.lines", c(list(x = densi$x, y = densi$y * contdensi), plot.line))
  } else {
    if(nb > 1)
      panel.rect(y = h$mids, x = 0, height = diff(breaks), width = heiBar,
                 col = plot.polygon$col, alpha = plot.polygon$alpha, border = plot.polygon$border, lty = plot.polygon$lty, 
                 lwd = plot.polygon$lwd, just = c("left", "center"), identifier = identifier)
    do.call("panel.lines", c(list(y = densi$x,  x = densi$y * contdensi), plot.line))
  }
}


adeg.panel.join <- function(drawLines, params = list()) {
  ## circle from c(0,0)p, radius = drawLines
  plot.line <- modifyList(list(add.line = trellis.par.get("add.line")), params, keep.null = TRUE)[[1L]] ## density line
  ## number of seg = 200
  plabels <- modifyList(adegpar("plabels"), params, keep.null = TRUE)[[1L]]
  scaleX <- c(current.panel.limits()$xlim, current.panel.limits()$ylim)
  xlines <- seq(from = min(scaleX) - 0.1 * min(scaleX), to = max(scaleX) * 1.1, length.out = 200)
  ylines <- lapply(drawLines, FUN = function(radius, x) {
    indx <- (x <= radius) ## x can be greated than radius
    return(c(sqrt(radius * radius - x[indx] * x[indx]), (- sqrt(abs(radius * radius - x[!indx] * x[!indx])))))
  }, x = xlines)
  
  trash <- lapply(ylines, FUN = function(y, x) {do.call("panel.lines", c(list(x = x[1:length(y)], y = y[1:length(y)]), plot.line))}, x = xlines)
  adeg.panel.label(x = sqrt(0.5) * drawLines, y = sqrt(0.5) * drawLines, as.character(drawLines), plabels)
}


## from http://rwiki.sciviews.org/doku.php?id=tips:graphics-grid:displaybitmap
## used in s.logo (rasterGrob) to handle pixmap objects

as.raster.pixmapRGB <- function(x, ...) {
  nr <- nrow(x@red)
  r <- rgb((x@red), (x@green), (x@blue))
  dim(r) <- x@size
  r
}


as.raster.pixmapGrey <- function(x, ...) {
  nr <- nrow(x@grey)
  r <- x@grey
  dim(r) <- x@size
  r
}
