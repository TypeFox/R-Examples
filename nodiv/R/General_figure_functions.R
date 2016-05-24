



##########################################################
# Here comes a list of functions to use for the analysis. These definitions should all be loaded into R


plot_grid <- function(x, coords, col, shape = NULL, shapefill = "grey", shapeborder = NA, zlim = NULL, zoom_to_points = FALSE, legend = TRUE, gridcol, gridlwd, gridsites, overlay_shape = FALSE, ...)
{
  if(inherits(x, "SpatialPixelsDataFrame"))
    rast <- raster(x) else
    if(missing(coords)) stop("coords must be defined if x is not a SpatialPixelsDataFrame") else
      if(inherits(coords, "SpatialPoints"))
      {
        if(!isGrid(coords)) stop("this function is only suitable for gridded data")
        coords <- SpatialPoints(coords)
        suppressWarnings(rast <- raster(SpatialPixelsDataFrame( coords, as.data.frame(x)))) #TODO NB
      } else {
        if(is.matrix(coords)) 
          coords <- as.data.frame(coords)
        if(is.data.frame(coords))
        {
          if(!is.vector(x)) stop("x must be a vector or a SpatialPixelsDataFrame")
          if(!nrow(coords) == length(x)) stop("x and coords must have the same number of elements")
          if(!ncol(coords) == 2) stop("coords must have exactly two columns, for the x and y coordinates of data")
          if(!isGrid(coords)) stop("this function is only suitable for gridded data")
          coords <- SpatialPoints(coords)
          suppressWarnings(rast <- raster(SpatialPixelsDataFrame(coords, as.data.frame(x)))) #TODO NB
        } else stop("Undefined arguments")
      }
  
  if(missing(col)){
    if(is.character(x)) x <- as.factor(x)
    if(is.factor(x))
      coltype <- "individual" else coltype <- "auto"      
    col <- choose.colors(getValues(rast), zlim, coltype = coltype)
    zlim <- attr(col, "zlim")
  }
  
  if(is.character(col) & length(col) == 1){
    if(col %in% c("auto", "ramp", "monochrome", "divergent", "individual")){
      val <- getValues(rast)
      if(col == "individual")
        val <- as.factor(as.character(val))
      col <- choose.colors(val, zlim, coltype = col)
      zlim <- attr(col, "zlim")
    }
  }
  
    if(is.character(col))
    if(length(col) == 1)
      if(!substr(col, 1, 1) == "#")
        if(!col %in% colors())
          col <- custom_palette(col)
  
  if(is.null(zlim)) zlim <- c(min(getValues(rast),na.rm = T),max(getValues(rast),na.rm = T))
  if(min(zlim) == 0 & identical(getValues(rast), floor(getValues(rast))) & length(col) > 1) col[1] <- "grey" #TODO an experimental hack
  

  
  if(is.null(shape)) plot(rast, zlim = zlim, col = col, legend = legend, ...) else
  {
    if(!inherits(shape, "SpatialPolygonsDataFrame"))
      shapeborder <- NULL
    if(zoom_to_points){
      if(inherits(shape, "Raster"))
        plot(shape, col = shapefill, border = shapeborder, xlim = bbox(coords)[1,], ylim = bbox(coords)[2,], legend = FALSE, ...) else 
          plot(shape, col = shapefill, border = shapeborder, xlim = bbox(coords)[1,], ylim = bbox(coords)[2,], ...)
    } else {
      if(inherits(shape, "Raster"))
        plot(shape, col = shapefill, border = shapeborder, legend = FALSE, ...) else plot(shape, col = shapefill, border = shapeborder, ...)
    }
    plot(rast, add = T, zlim = zlim, col = col, legend = legend)
  }
  
  if(!missing(gridcol) || !missing(gridlwd) || !missing(gridsites)){
    if(missing(gridcol))
      gridcol <- "lightgrey"
    if(missing(gridlwd))
      gridlwd <- 0.4
    if(missing(gridsites))
      gridsites <- 1:length(x)
    add_grid(coords, sites = gridsites, border = gridcol, lwd = gridlwd)
  }
  
  if(overlay_shape){
    if(is.null(shape))
      stop("No shape to overlay in the data set")
    if(!inherits(shape, "SpatialPolygonsDataFrame"))
      stop("Can only overlay objects of type SpatialPolygonsDataFrame")
    if(is.na(shapeborder))
      shapeborder = "white"
    plot(shape, border = shapeborder, add = TRUE, lwd = 0.15, ...) 
  }
    
  
  invisible(rast)
}


plot_points <- function(x, coords, col , shape = NULL, shapefill = "grey", zlim= NULL,  zoom_to_points = FALSE, pch = 16, bg = par("bg"), legend = TRUE, ...)
{  
  if(inherits(x, "SpatialPointsDataFrame"))
  {
    coords <- SpatialPoints(x)
    x <- as.numeric(x@data[,1])
  } else {
      if(missing(coords)) stop("coords must be defined if x is not a SpatialPixelsDataFrame") else {
        if(!inherits(coords, "SpatialPoints"))
        {
          if(is.matrix(coords)) coords <- as.data.frame(coords) 
          if(is.data.frame(coords))
          {
            if(!is.vector(x)) stop("x must be a vector or a SpatialPointsDataFrame")
            if(!nrow(coords) == length(x)) stop("x and coords must have the same number of elements")
            if(!ncol(coords) == 2) stop("coords must have exactly two columns, for the x and y coordinates of data")
            coords <- SpatialPoints(coords)
            
          } else stop("Undefined arguments")
        }
    }
  }

  if(missing(col)){
    col <- choose.colors(x, zlim)
    zlim <- attr(col, "zlim")
  }
  
  if(is.null(zlim)) zlim <- c(min(x,na.rm = T),max(x,na.rm = T))
  
  
  coords <- SpatialPoints(coords)

  #split.screen( rbind(c(0, .8,0,1), c(.8,1,0,1)))
  #screen(1)
  
  oldpar <- par()
  #par(plt = c(0,0.8,0,1))
  par(mar = c(5,4,4,6) + 0.1)
  
  if(is.character(col))
    if(length(col) == 1)
      if(!substr(col, 1, 1) == "#")
        if(!col %in% colors())
          col <- custom_palette(col)
    
  plotcol <- create.cols(x, col, zlim = zlim)
  if(pch %in% 21:25)
  {
    bg <- plotcol
    plotcol <- 1
  }
  coords <- coords[order(x),]
  plotcol <- plotcol[order(x)]  
  
  if(is.null(shape)) plot(coords, col = plotcol, pch = pch, bg = bg, ...) else
  {
    if(inherits(shape, "Raster")) legend <- FALSE else legend <- NULL
    if(inherits(shape, "SpatialPolygonsDataFrame"))
      border <- shapefill else border <- NULL
    if(zoom_to_points)
      plot(shape, col = shapefill,  xlim = bbox(coords)[1,], ylim = bbox(coords)[2,], legend = legend, border = border, ...) else  plot(shape, col = shapefill, legend = legend, border = border, ...)
    plot(coords, col = plotcol, pch = pch, bg = bg, add = T)
  } 

  #screen(2)
  par <- oldpar


  if(legend)  add_legend(col = col, zlim = zlim)
  #image.plot( zlim = zlim,legend.only=TRUE, smallplot=c(.85,.87, .38,.65), col=col)
  ret <- data.frame(x)
  names(ret) <- deparse(substitute(x))
  invisible(SpatialPointsDataFrame(coords, ret))
}


plot_nodes_phylo <- function(variable, tree, label = variable, main = deparse(substitute(variable)), zlim = NULL, col , show.legend = TRUE, sig.cutoff, nodes, roundoff= TRUE, show.tip.label = NULL, cex = NULL, ...)
{
    if(!length(variable) == Nnode(tree))
    stop("The length of the variable vector must be the same length as the number of nodes on the tree")
  
  if(is.null(show.tip.label)) show.tip.label <- isTRUE(Ntip(tree) < 50)
  
  plotvar <- variable
  if(roundoff & is.numeric(label)) label = round(label,2)
  

  if(missing(col)){
    col <- choose.colors(plotvar, zlim)
    zlim <- attr(col, "zlim")
  }
  
  if(is.character(col))
    if(length(col) == 1)
      if(!substr(col, 1, 1) == "#")
        if(!col %in% colors())
          col <- custom_palette(col)
  
  
  if(is.null(zlim)) zlim <- c(min(plotvar, na.rm = T), max(plotvar, na.rm = T))
  
  if(is.null(cex)) cex <- par("cex")
  sizes <- cex * 4 * sqrt((plotvar - zlim[1])/zlim[2])
  
  if(missing(nodes)) node_index = rep(TRUE, Nnode(tree)) else {
    node_index <- rep(FALSE, Nnode(tree))
    node_index[nodes] <- TRUE
  }
  
  node_index[is.na(variable)] <- FALSE
  
  if(!missing(sig.cutoff))
    node_index[variable <= sig.cutoff] <- FALSE
  
  nodes <- nodenumbers(tree)[node_index]
  plotvar <- plotvar[node_index]
  label <- label[node_index]
  sizes <- sizes[node_index]
  
  
  if(show.legend)
  {
  #  split.screen( rbind(c(0, .9,0,1), c(.9,1,0,1)))
  #  screen(1)
    oldpar <- par()
    par(mar = c(5,4,4,6) + 0.1)
  }
  
  plot(tree, show.tip.label = show.tip.label, ...)
  title(main)
  nodelabels(pch = 16, node = nodes, col = create.cols(plotvar, col, zlim), cex = sizes)
  nodelabels(text = as.character(label), node = nodes, cex = 0.6, frame = "none") 
  
  if(show.legend)
  {
    par <- oldpar #screen(2)
    add_legend(col = col, zlim = zlim)
  }
}

add_legend <- function (zlim, smallplot=c(.85,.866, .38,.65), col)
{
  old.par <- par()
  if ((smallplot[2] < smallplot[1]) | (smallplot[4] < smallplot[3])) {
    stop("plot region too small to add legend\n")
  }
  
  ix <- 1
  minz <- zlim[1]
  maxz <- zlim[2]
  nlevel <- length(col)
  binwidth <- (maxz - minz)/nlevel
  midpoints <- seq(minz + binwidth/2, maxz - binwidth/2, by = binwidth)
  iy <- midpoints
  iz <- matrix(iy, nrow = 1, ncol = length(iy))
  
  par(new = TRUE, pty = "m", plt = smallplot, err = -1)
  
  axis.args <- list(side =  4, mgp = c(3, 1, 0), las = 2)
  image(ix, iy, iz, xaxt = "n", yaxt = "n", xlab = "", ylab = "", col = col)
  
  do.call("axis", axis.args)
  box()
  par(new = FALSE, pty = old.par$pty, plt = old.par$plt, err = old.par$err)
  invisible()
}

add_grid <- function(coords, sites = 1:nrow(coords), border = "lightgrey", lwd = 0.4, ...){
  if(inherits(coords, "SpatialPoints"))
    coords <- coordinates(coords)
  lon <- coords[,1]
  lat <- coords[,2]	
  cellsize_x <- min(abs(diff(unique(lon))))/2
  cellsize_y <- min(abs(diff(unique(lat))))/2
  lon <- lon[sites]
  lat <- lat[sites]
  rect(lon - cellsize_x, lat - cellsize_y, lon + cellsize_x, lat + cellsize_y, col = NA, border = border, lwd = lwd, ...)
}

phyplot <- function(x, type = "phylogram", use.edge.length = TRUE, node.pos = NULL, 
          show.tip.label = TRUE, show.node.label = FALSE, edge.color = "black", 
          edge.width = 1, edge.lty = 1, font = 3, cex = par("cex"), 
          adj = NULL, srt = 0, no.margin = FALSE, root.edge = FALSE, 
          label.offset = 0, underscore = FALSE, x.lim = NULL, y.lim = NULL, 
          direction = "rightwards", lab4ut = NULL, tip.color = "black", 
          plot = TRUE, rotate.tree = 0, open.angle = 0, node.depth = 1, 
          align.tip.label = FALSE, ...) 
  plot(x, type, use.edge.length, node.pos, 
       show.tip.label, show.node.label, edge.color, 
       edge.width, edge.lty, font, cex, 
       adj, srt, no.margin, root.edge, 
       label.offset, underscore, x.lim, y.lim, 
       direction, lab4ut, tip.color, 
       plot, rotate.tree, open.angle, node.depth, 
       align.tip.label, ...) 

