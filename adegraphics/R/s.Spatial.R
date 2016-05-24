s.Spatial <- function(spObj, col = TRUE, nclass = 5, plot = TRUE, storeData = TRUE, pos = -1, ...) {
  oldparamadeg <- adegpar()
  on.exit(adegpar(oldparamadeg))
  sortparameters <- sortparamADEg(...)
  adegtot <- adegpar(sortparameters$adepar)
  
  xy.spObj <- coordinates(spObj)[, , drop = FALSE]  ## to access 'coordinates' in the 'imports' environment of 'adegraphics'
  
  ## default values for non-used parameters
  defaultpar <- list(plabels = list(cex = 0), pgrid = list(draw = FALSE), ppoints = list(cex = 0), porigin = list(include = FALSE))
  sortparameters$adepar <- modifyList(defaultpar, sortparameters$adepar, keep.null = TRUE)
  if(is.logical(col)){
      if(col)
          colnew <- adegtot$pSp$col
      else
          colnew <- "transparent"	## col == FALSE
  } else{
      colnew <- col
  }
  
  nvar <- 0
  if(length(grep("DataFrame", class(spObj))) > 0)
    nvar <- ncol(spObj)
  
  ## limits management 
  limsSp <- bbox(spObj)
  lim.global <- setlimits2D(minX = limsSp[1, 1], maxX = limsSp[1, 2], minY = limsSp[2, 1], maxY = limsSp[2, 2], includeOr = FALSE) 
	if(is.null(sortparameters$g.args$xlim))
  	sortparameters$g.args$xlim <- lim.global$xlim
  if(is.null(sortparameters$g.args$ylim))
    sortparameters$g.args$ylim <- lim.global$ylim
  
  if(nvar < 2) {
    if(nvar == 1) {
      ## Spatial*DataFrame object -> ADEg
      defaultpar <- list(psub = list(text = names(spObj)[1]))
      sortparameters$adepar <- modifyList(defaultpar, sortparameters$adepar, keep.null = TRUE)
      if(is.logical(col)) {
        if(col) {
          if(is.numeric(spObj@data[, 1])) {
            nclasspretty <- length(pretty(spObj@data[, 1], nclass)) - 1
            nclasspretty <- length(pretty(spObj@data[, 1], nclasspretty)) - 1 ## repeated in order to have always the same number of class
            if(is.null(sortparameters$adepar$pSp$col))
              colnew <- adegtot$ppalette$quanti(nclasspretty)
          } else
            if(is.null(sortparameters$adepar$pSp$col))
              colnew <- adegtot$ppalette$quali(nlevels(as.factor(spObj@data[, 1])))
        }
      } 
    } 
    
    sortparameters$adepar$pSp$col <- colnew

    ## create map
    object <- do.call("s.label", c(list(dfxy = xy.spObj, Sp = substitute(spObj), plot = FALSE, storeData = storeData, pos = pos - 2), sortparameters$adepar, sortparameters$trellis, sortparameters$g.args))
    
  } else {
    ## Spatial*DataFrame object with several variables -> ADEgS
    listGraph <- list()
    for(i in 1:nvar) {
      sortparameters$adepar <- modifyList(defaultpar, sortparameters$adepar, keep.null = TRUE)
      if(is.logical(col)) {
        if(col) {
          if(is.numeric(spObj@data[, i])) {
            nclasspretty <- length(pretty(spObj@data[, i], nclass)) - 1
            nclasspretty <- length(pretty(spObj@data[, i], nclasspretty)) - 1 ## repeated in order to have always the same number of class
            if(is.null(sortparameters$adepar$pSp$col))
              colnew <- adegtot$ppalette$quanti(nclasspretty)
          } else
            if(is.null(sortparameters$adepar$pSp$col))
              colnew <- adegtot$ppalette$quali(nlevels(as.factor(spObj@data[, i])))
        }
      } else {
        colnew <- col
      }
      
      sortparameters$adepar$pSp$col <- colnew
      sortparameters$adepar$psub$text <- names(spObj)[i]
      
      ## create map
      listGraph <- c(listGraph, do.call("s.label", c(list(dfxy = xy.spObj, Sp = substitute(spObj[, i]), plot = FALSE, storeData = storeData, pos = pos - 2), sortparameters$adepar, sortparameters$trellis, sortparameters$g.args)))
    }
    names(listGraph) <- names(spObj)
    posmatrix <- layout2position(.n2mfrow(nvar), ng = nvar)
    object <- new(Class = "ADEgS", ADEglist = listGraph, positions = posmatrix, add = matrix(0, ncol = nvar, nrow = nvar), Call = match.call())
  }
  
  if(plot)
    print(object)
  invisible(object) 
}

