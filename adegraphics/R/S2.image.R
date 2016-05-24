##########################################################################
##                           s.image                                    ##
##########################################################################

## TODO: prendre en comptre les differents z
setClass(
  Class = "S2.image",
  contains = "ADEg.S2",
)


setMethod(
  f = "initialize",
  signature = "S2.image",
  definition = function(.Object, data = list(dfxy = NULL, z = NULL, xax = 1, yax = 2, frame = 0, storeData = TRUE), ...) {
    .Object <- callNextMethod(.Object, data = data, ...) ## ADEg.S2 initialize
    .Object@data$z <- data$z
    return(.Object)
  })


setMethod(
  f = "prepare",
  signature = "S2.image",
  definition = function(object) {
    name_obj <- deparse(substitute(object))
    
    ## pre-management of graphics parameters      
    oldparamadeg <- adegpar()
    on.exit(adegpar(oldparamadeg))
    adegtot <- adegpar(object@adeg.par)
   
    object@g.args$gridsize <- rep(object@g.args$gridsize, length.out = 2)
    if(is.null(object@adeg.par$porigin$include) & (any(names(object@g.args) %in% c("Sp", "nbobject", "outsideLimits"))))
      adegtot$porigin$include <- FALSE

    if(object@data$storeData) {
      dfxy <- object@data$dfxy
      z <- object@data$z
    } else {
      dfxy <- eval(object@data$dfxy, envir = sys.frame(object@data$frame))
      z <- eval(object@data$z, envir = sys.frame(object@data$frame))
    }
    
    if(is.null(object@g.args$breaks))
        object@s.misc$breaks.update <- pretty(z, object@g.args$nclass)
    else
        object@s.misc$breaks.update <- object@g.args$breaks
    
    object@s.misc$breaks.update <- breakstest(object@s.misc$breaks.update, z, n = length(object@s.misc$breaks.update))
    n <- length(object@s.misc$breaks.update)
    
    if(!is.null(object@g.args$col)) {
        if(length(object@g.args$col) < (n - 1))
            stop(paste("not enough colors defined, at least ", (n - 1), " colors expected", sep = ""), call. = FALSE)
        adegtot$ppoints$col <- object@g.args$col[1:(n - 1)]  ## color given by the user
    } else {
        if(is.null(object@adeg.par$ppoints$col))
            adegtot$ppoints$col <- adegtot$ppalette$quanti(n - 1)
    }

    ## object modification before calling inherited method
    object@adeg.par <- adegtot
    callNextMethod() ## prepare graph (provide limits that are used below)
    
    ## create a sp grid
    minX <- object@g.args$xlim[1]
    minY <- object@g.args$ylim[1]
    cgridX <- diff(object@g.args$xlim) / object@g.args$gridsize[1]
    cgridY <- diff(object@g.args$ylim) / object@g.args$gridsize[2]
    gridSp <- SpatialGrid(GridTopology(c(minX, minY), c(cgridX, cgridY), c(object@g.args$gridsize[1], object@g.args$gridsize[2])))

    x <- dfxy[, object@data$xax]
    y <- dfxy[, object@data$yax]
    if(!is.null(object@g.args$outsideLimits)) {
      ## outside limits are provided by an sp object
      whichis <- over(gridSp, object@g.args$outsideLimits)
    } else {
      ## define the outside limits by convex hull
      beplot <- cbind(x, y)[chull(cbind(x, y)), ] 
      extCoord <- SpatialPolygons(list(Polygons(list(Polygon(rbind(cbind(beplot[, 1], beplot[, 2]), beplot[1, ]))), ID = "extcoord")))
      whichis <- over(gridSp, extCoord)
    }
    
    ## NA not handled by panel.levelplot call afterward ==> we remove the points
    newgrid <- coordinates(gridSp)
    names(newgrid) <- c("x", "y")
    lo <- loess(z ~ x + y, span = object@g.args$span) ## Local Polynomial Regression Fitting
    predictval <- predict(lo, newdata = newgrid)            
    predictval[which(is.na(whichis))] <- NA
    tokeep <- !is.na(predictval)
    predictval <- predictval[tokeep]
    newgrid <- newgrid[tokeep, ]
    object@stats$value <- predictval
    object@s.misc$newgrid <- newgrid
      
    assign(name_obj, object, envir = parent.frame())
  })


setMethod(
  f = "panel",
  signature = "S2.image",
  definition = function(object, x, y) {
    zvalue <- object@stats$value
    col <- object@adeg.par$ppoints$col
    xx <- object@s.misc$newgrid[, 1]
    yy <- object@s.misc$newgrid[, 2]
    panel.levelplot(x = xx, y = yy, z = zvalue, subscripts = TRUE, col.regions = col, contour = object@g.args$contour, region = object@g.args$region, labels = object@adeg.par$plabels,
                    label.style = if(object@adeg.par$plabels$srt == "horizontal") "flat" else "align")
  })


s.image <- function(dfxy, z, xax = 1, yax = 2, span = 0.5, gridsize = c(80L, 80L), contour = TRUE, region = TRUE, outsideLimits = NULL, breaks = NULL, nclass = 8, 
  									col = NULL, facets = NULL, plot = TRUE, storeData = TRUE, add = FALSE, pos = -1, ...) {
 
  if(!is.null(outsideLimits)) {
    if(!inherits(outsideLimits, "SpatialPolygons"))
      stop("limits must be a SpatialPolygons")
  }
  
  ## evaluation of some parameters (required for multiplot)
  thecall <- .expand.call(match.call())
  df <- try(as.data.frame(eval(thecall$dfxy, envir = sys.frame(sys.nframe() + pos))), silent = TRUE)
  z <- eval(thecall$z, envir = sys.frame(sys.nframe() + pos))
  if((class(df) == "try-error") | is.null(thecall$dfxy)) ## non convenient dfxy argument
    stop("non convenient selection for dfxy (can not be converted to dataframe)")
  if(NROW(df) != NROW(z))
    stop("dfxy and z should have the same number of rows")

  ## parameters sorted
  sortparameters <- sortparamADEg(...)
  
  ## facets
  if(!is.null(facets)) {
    if((length(xax) == 1 & length(yax) == 1) & NCOL(z) == 1)
      object <- multi.facets.S2(thecall, sortparameters$adepar, samelimits = sortparameters$g.args$samelimits)
    else 
      stop("Facets are not allowed with multiple xax/yax or multiple z")
  }
  
  ## multiple axes
  else if((length(xax) > 1 | length(yax) > 1)) {
    if(NCOL(z) == 1)
      object <- multi.ax.S2(thecall)
    else 
      stop("Multiple xax/yax are not allowed with multiple z")
  }
  
  ## multiple z
  else if(NCOL(z) > 1) {
    object <- multi.variables.S2(thecall, "z")
  }
  
  ## simple ADEg graphic
  else {
    if(length(sortparameters$rest))
      warning(c("Unused parameters: ", paste(unique(names(sortparameters$rest)), " ", sep = "")), call. = FALSE)
    
    ## creation of the ADEg object
    g.args <- c(sortparameters$g.args, list(breaks = breaks, nclass = nclass, span = span, gridsize = gridsize, outsideLimits = outsideLimits,
                                            contour = contour, region = region, col = col))
    if(storeData)
    	tmp_data <- list(dfxy = dfxy, xax = xax, yax = yax, z = z, frame = sys.nframe() + pos, storeData = storeData)
    else
    	tmp_data <- list(dfxy = thecall$dfxy, xax = xax, yax = yax, z = thecall$z, frame = sys.nframe() + pos, storeData = storeData)
    object <- new(Class = "S2.image", data = tmp_data, adeg.par = sortparameters$adepar, trellis.par = sortparameters$trellis, g.args = g.args, Call = as.call(thecall))
   
    ## preparation of the graph
    prepare(object)
    setlatticecall(object)
    if(add)
      object <- add.ADEg(object)
  }
  
  if(! add & plot)
    print(object)
  invisible(object)
}

