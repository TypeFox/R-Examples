#' @import grid
NULL

#' sliceData
#' 
#' extract a 2D slice from a \code{BrainVolume} instance.
#' 
#' @param vol an \code{BrainVolume} instance
#' @param slice the integer index of the slice to cut.
#' @param axis the axis number (1, 2, 3) defining fixed axis of the 2D slice.
#' @export
sliceData <- function(vol, slice, axis=3) {
  imslice <- switch(axis,
                    "1"=vol[slice,,],
                    "2"=vol[,slice,],
                    "3"=vol[,,slice])
  
  imslice <- t(imslice[nrow(imslice):1, ncol(imslice):1,drop=FALSE])    
  
  
}

#' mapToColors
#' 
#' map an matrix of intensity values to a matrix of color values.
#' 
#' @importFrom grDevices heat.colors
#' @param imslice an image matrix defining intensity values
#' @param col a color map
#' @param zero.col the background color.
#' @param alpha transparency multiplier
#' @importFrom grDevices col2rgb
#' @export
mapToColors <- function(imslice, col=heat.colors(128, alpha = 1), zero.col = "#00000000", alpha=1) {
  vrange <- range(imslice)
  imcols <- col[(imslice - vrange[1])/diff(vrange) * (length(col) -1) + 1]
  
  dim(imcols) <- dim(imslice)
  imcols[imslice == 0] <- zero.col
  
  if (alpha < 1) {
    rgbmat <- col2rgb(imcols, alpha=TRUE)
    rgbmat <- rgbmat/255
  
    if (alpha < 1) {
      rgbmat[4,] <- rgbmat[4,] * alpha
    }
  
    array(t(rgbmat), c(dim(imslice), 4))
  } else {
    imcols
  }
  
}

#' image
#' @param slice the voxel index of the slice to display
#' @param col a color map
#' @param zero.col the color to use when the value is 0 (e.g background color)
#' @param ... extra arguments to passed to \code{grid.raster}
#' @export
#' @rdname image-methods
setMethod(f="image", signature=signature(x = "BrainVolume"),
          def=function(x, slice, col=gray((0:255)/255, alpha=1), zero.col = "#000000", axis=3, ...) {    
            imslice <- sliceData(x, slice, axis)
            imcols <- mapToColors(imslice, col, zero.col)
            ras <- as.raster(imcols)
            ras[imslice == 0] <- zero.col
            
            grid.newpage()
            grid.raster(ras, ...)
          })

#' Layer
#' 
#' create a \code{\linkS4class{Layer}} object
#' @param vol volume instance of \code{\linkS4class{BrainVolume}}
#' @param colorMap a lookup table defining mapping from image intensity values to colors.
#' @param thresh a range (min,max) defining the threshold window for determining image opacity.
#' @param axis the axis index of the axis perpendicular to the xy plane (options: 1,2,3; default is 3)
#' @param zero.col the color used when the value is zero.
#' @param alpha transparency multiplier, vlaue between 0 and 1.
#' @return an object of class \code{Layer}
#' @export
#' @rdname Layer
#' @importFrom grDevices gray
Layer <- function(vol, colorMap=gray((0:255)/255, alpha=1), thresh=c(0,0), axis=3, zero.col="#000000", alpha=1) {
  new("Layer", vol=vol, colorMap=colorMap, thresh=thresh, axis=axis, zero.col=zero.col, alpha=alpha)
}




#' as.raster
#' 
#' @export 
#' @param x the layer to convert
#' @param zpos the z coordinate in coordinate space
#' @rdname as.raster-methods
setMethod(f="as.raster", signature=signature(x = "Layer"),
          def=function(x, zpos) {  
            slice <- axisToIndex(space(x@vol), zpos, x@axis)
            imslice <- sliceData(x@vol, slice, x@axis)     
            vrange <- range(imslice)
            
            thresh <- x@thresh
            
            lookup <- (imslice - vrange[1])/diff(vrange) * (length(x@colorMap) -1) + 1
            imcols <- x@colorMap[lookup]
            
            if (length(thresh) == 1) {
              thresh <- c(-Inf, thresh)
            }
                      
            imcols[imslice == 0] <- x@zero.col
            
            if (diff(thresh) > 0) {
              imcols[(imslice >= thresh[1] & imslice <= thresh[2])] <- "#00000000"
            }
            
            dim(imcols) <- dim(imslice)
            ras <- as.raster(imcols)
            #ras[imslice == 0] <- zero.col            
            ras
          })




#' overlay
#' 
#' @export 
#' @rdname overlay-methods
setMethod(f="overlay", signature=signature(x = "Layer", y="Layer"),
          def=function(x, y) {  
            new("Overlay", layers=list(x,y))  
          })


#' @export 
#' @rdname overlay-methods
#' @param e1 the left operand
#' @param e2 the right operand
setMethod(f="+", signature=signature(e1 = "Overlay", e2="Layer"),
          def=function(e1, e2) {  
            new("Overlay", layers=c(e1@layers, e2))
          })

#' @export 
#' @rdname overlay-methods
setMethod(f="+", signature=signature(e1 = "Layer", e2="Layer"),
          def=function(e1, e2) {  
            new("Overlay", layers=list(e1, e2))
          })

#' image
#' @param x the object to display
#' @param zpos the z coordinate
#' @param axis the axis index
#' @export
#' @rdname image-methods
setMethod(f="image", signature=signature(x = "Overlay"),
          def=function(x, zpos, axis=3) {  
            #grid.newpage()
            for (layer in x@layers) {
              ras <- as.raster(layer, zpos, layer@thresh, axis=axis) 
              grid.raster(ras, interpolate=TRUE)
            }             
          })


#' @rdname image-methods
#' @export
setMethod(f="image", signature=signature(x = "Layer"),
          def=function(x, zpos, axis=3) {  
            ras <- as.raster(x, zpos, x@thresh,axis=axis)      
            #grid.newpage()
            grid.raster(ras, interpolate=TRUE)
          })


#' @export
#' @rdname renderSlice-methods
#' @param zero.col color used when background intensity is 0.
#' @param units grid unit type, e.g. "mm", "inches"
setMethod(f="renderSlice", signature=signature(x="Overlay", zpos="numeric", width="numeric", height="numeric", colmap="missing"),
          def=function(x, zpos, width, height, zero.col="#000000FF", units="mm") {
            sliceList <- lapply(x@layers, function(layer) {
              renderSlice(layer, zpos=zpos, width=width, height=height, zero.col=zero.col, units=units)
            })
            
            slices <- lapply(sliceList, function(x) x@slice)
            grobs <- lapply(sliceList, function(x) x@raster)
            gl <- do.call(gList, grobs)
            
            new("RenderedSliceStack", slices=slices, width=width, height=height, grob=gl)
            
          })


#' @export
#' @rdname renderSlice-methods
setMethod(f="renderSlice", signature=signature(x="Layer", zpos="numeric", width="numeric", height="numeric", colmap="missing"),
          def=function(x, zpos, width, height, colmap, zero.col="#000000FF", units="mm") {
            slice <- slice(x@vol, axisToIndex(space(x@vol), zpos, x@axis), x@axis, "")
            grob <- render(slice, width, height, colmap=x@colorMap, zero.col=zero.col, alpha=x@alpha, units=units)
            
            new("RenderedSlice", slice=slice, width=width, height=height, raster=grob)
          })

#' @export
#' @rdname render-methods
#' @param zero.col color used when background intensity is 0.
#' @param alpha transparency multiplier
#' @param units grid unit type, e.g. "mm", "inches"
setMethod(f="render", signature=signature(x="BrainSlice", width="numeric", height="numeric", colmap="character"),
          def=function(x, width, height, colmap, zero.col="#000000FF", alpha=1, units="mm") {
            imslice <- t(x@.Data[nrow(x@.Data):1, ncol(x@.Data):1,drop=FALSE])    
            imcols <- mapToColors(imslice, colmap, zero.col, alpha=alpha)
            ras <- as.raster(imcols)
  
            grob <- rasterGrob(ras, width=unit(width, units), height=unit(height, units), interpolate=TRUE)
    
          })



#orthoPlot <- function(layer, zpos) {
#  
#  vpmain <- viewport(x=0, y=0, width=1, height=1)
#  vptopright <- viewport(x=unit(.5, "npc"), y=unit(0, "npc"), width=unit(.5, "npc"), height=unit(.5, "npc"))
#  vpbottomright <- viewport(x=unit(.5, "npc"), y=unit(.5, "npc"), width=unit(.5, "npc"), height=unit(.5, "npc"))
#  vpleft <- viewport(x=unit(0, "npc"), y=unit(0,"ncp"), width=unit(.5, "npc"), height=unit(1, "npc"))
#}

# plotMontage <- function(x, layout=c(3,3), zstart, zend) {
#   
#   zslices <- seq(zstart, zend,by=(zend-zstart)/prod(layout))
#   
#   raslist <- lapply(zslices, function(z) {
#     print(z)
#     zind <- axisToIndex(space(x@vol), z, x@axis)
#     dat <- slice(x@vol, zind, x@axis,"")
#     xy <- indexToCoord(space(dat), 1:length(dat))
#     list(xy=as.matrix(xy), values=as.numeric(dat@.Data), slice=z)
#   })
#   
#   xy <- do.call(rbind, lapply(raslist, "[[", "xy"))
#   dfras <- data.frame(x=xy[,1], y=xy[,2], values=unlist(lapply(raslist, "[[", "values")), slice=sapply(raslist, "[[", "slice"))
#   dfras$slice <- factor(dfras$slice)
#   
#   p <- ggplot(data=dfras, aes(x=x,y=y)) +
#     theme_bw() + coord_equal() +
#     geom_raster(aes(fill=values)) +
#     facet_grid(. ~ slice)
#   
#     
#   
# }
# 

#' imageGrid
#' 
#' Display a set of images slices in a 2D montage
#' 
#' @param layer the layer to display
#' @param gridDim the dimensions of the 2D grid montage
#' @param zstart the z coordinate of the first slice
#' @param zend the z coordinate of the last slice
#' @param panelSize the size of each panel in the montage (default unit is inches)
#' @param panelUnit the unit for the panel size (default is inches)
#' @param interpolate whether to interpolate pixel values
#' @param fontCol color of labels indicating slice level
#' @rdname imageGrid
imageGrid <- function(layer, gridDim=c(3,3), zstart, zend, panelSize=3, panelUnit="inches", interpolate=FALSE, fontCol="red") {
  slices <- seq(zstart, zend, length.out=prod(gridDim))
  grid.newpage()
  layout <- grid.layout(gridDim[1], gridDim[2], widths=rep(unit(panelSize, panelUnit), gridDim[2]),  
                                                heights=rep(unit(panelSize, panelUnit), gridDim[1]))
  
  grid.rect(unit(0, "npc"), y=unit(0, "npc"), just=c("left", "bottom"), gp=gpar(fill="black"))
  pushViewport(viewport(layout=layout))
  
  scount = 1
  for (i in 1:gridDim[1]) {
    for (j in 1:gridDim[2]) {
      pushViewport(viewport(layout.pos.col=j, layout.pos.row=i))
      ras <- as.raster(layer, slices[scount])
      
      grid.raster(ras, interpolate=interpolate)
      grid.text(paste(round(slices[scount])), x=unit(.15, "npc"), y=unit(.1, "npc"), 
                just="centre", gp=gpar(fontsize=14, col=fontCol))
      popViewport()
      scount <- scount+1
    }
      
  }
  
}
          

