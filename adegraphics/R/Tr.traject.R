######################################################
##                   Tr.traject                    ###
######################################################

setClass(
  Class = "Tr.traject",
  contains = "ADEg.Tr"
)


setMethod(
  f = "initialize",
  signature = "Tr.traject",
  definition = function(.Object, data = list(dfxyz = NULL, fac = NULL, labels = NULL, frame = 0, storeData = TRUE), ...) {
    .Object <- callNextMethod(.Object, data = data, ...)
    .Object@data$fac <- data$fac
    .Object@data$labels <- data$labels
    return(.Object)
  })


setMethod(
  f = "prepare",
  signature = "Tr.traject",
  definition = function(object) {
    name_obj <- deparse(substitute(object))
    
    if(object@data$storeData) {
      df <- object@data$dfxyz
      fac <- as.factor(object@data$fac)
    } else {
      df <- eval(object@data$dfxyz, envir = sys.frame(object@data$frame))
      fac <- as.factor(eval(object@data$fac, envir = sys.frame(object@data$frame)))
    }
    
    ## pre-management of graphics parameters      
    oldparamadeg <- adegpar()
    on.exit(adegpar(oldparamadeg))
    adegtot <- adegpar(object@adeg.par)
    
    ## change default for some parameters
    if(!is.null(object@g.args$col))
      if(is.logical(object@g.args$col)) {
        if(object@g.args$col)
          adegtot$ppoints$col <- adegtot$ppoints$fill <- adegtot$plabels$col <- adegtot$plabels$boxes$border <- adegtot$plines$col <- adegtot$ppalette$quali(nlevels(fac))
      }  
      else
        adegtot$ppoints$col <- adegtot$ppoints$fill <- adegtot$plabels$col <- adegtot$plabels$boxes$border <- adegtot$plines$col <- rep(object@g.args$col, length.out = nlevels(fac))
    
    ## object modification before calling inherited method
    object@adeg.par <- adegtot
    callNextMethod() ## prepare graph
    
    ## calculate 2D coordinates
    df <- sweep(df, 1, rowSums(df), "/")
    object@stats$coords2d <- .coordtotriangleM(df, mini3 = object@g.args$min3d, maxi3 = object@g.args$max3d)[, 2:3]

    object@adeg.par$plabels$optim <- FALSE
    assign(name_obj, object, envir = parent.frame())
  })


setMethod(
  f = "panel",
  signature = "Tr.traject",
  definition = function(object, x, y) {
    
    if(object@data$storeData) {
      fact <- object@data$fac
      labels <- object@data$labels
    } else {
      fact <- eval(object@data$fac, envir = sys.frame(object@data$frame))
      labels <- eval(object@data$labels, envir = sys.frame(object@data$frame))
    }
    
    todrawX <- split(object@stats$coords2d[, 1], fact)
    todrawY <- split(object@stats$coords2d[, 2], fact)
    sizelevels <- sapply(todrawX, length)
    if(!is.null(object@g.args$order))
      orderdraw <- split(order, fact)
    else
      orderdraw <- lapply(sizelevels, FUN = function(x) if(x > 0)  1:x else NULL)
    ## ordrerdraw is a list used to recycle graphical parameters
    
    setparam <- function(params, nblevel, sizelevels) {
      ## for param begin and end or repetition                               
      if(length(params) == nblevel)
        return(mapply(params, FUN = function(x, y) rep(x, length.out = y), sizelevels, SIMPLIFY = FALSE))
      else
        return(mapply(sizelevels, FUN = function(x, y) rep(params, length.out = x), SIMPLIFY = FALSE))    
    }
       
    parrows <- lapply(object@adeg.par$parrows, setparam, nblevel = length(todrawX), sizelevels = sizelevels)
    plines <- lapply(object@adeg.par$plines, setparam, nblevel = length(todrawX), sizelevels = sizelevels)
    ppoints <- lapply(object@adeg.par$ppoints, setparam, nblevel = length(todrawX), sizelevels = sizelevels)
    
    for(i in 1:length(todrawX)) {
      if(length(todrawX[[i]]) > 0)
        panel.points(x = todrawX[[i]], y = todrawY[[i]], col = ppoints$col[[i]], cex = ppoints$cex[[i]], pch = ppoints$pch[[i]], fill = ppoints$fill[[i]])
    }
    
    for(i in 1:length(todrawX)) {
      if(length(todrawX[[i]]) > 1) {
        suborder <- orderdraw[[i]]
        for(j in 1:(length(todrawX[[i]]) - 1)) {
          panel.arrows(x0 = todrawX[[i]][suborder[j]], y0 = todrawY[[i]][suborder[j]],
            x1 = todrawX[[i]][suborder[j + 1]], y1 = todrawY[[i]][suborder[j + 1]],
            angle = parrows$angle[[i]][suborder[j + 1]], length = parrows$length[[i]][suborder[j + 1]],
            ends = parrows$end[[i]][suborder[j + 1]], lwd = plines$lwd[[i]][suborder[j + 1]],
            col = plines$col[[i]][suborder[j + 1]], lty = plines$lty[[i]][suborder[j + 1]])
        }
      }
    }
    
    if(any(object@adeg.par$plabels$cex > 0)) {
      ## draws labels in the middle part of the trajectory
      middl <- sapply(orderdraw, FUN = function(x) floor(length(x) / 2))
      x <- y <- rep(NA, length(middl))
      for(i in 1:length(middl)) {
        if(length(todrawX[[i]]) > 1) {
          x[i] <- (todrawX[[i]][suborder[middl[i]]] + todrawX[[i]][suborder[middl[i]+1]]) / 2
          y[i] <- (todrawY[[i]][suborder[middl[i]]] + todrawY[[i]][suborder[middl[i]+1]]) / 2
        }
      }
      
      ## always false here no optimization for straject object
      object@adeg.par$plabels$optim <- FALSE
      adeg.panel.label(x, y, labels = labels, plabels = object@adeg.par$plabels)
    }
  })


triangle.traject <- function(dfxyz, fac = gl(1, nrow(dfxyz)), order, labels = levels(fac), col = NULL, adjust = TRUE, 
  min3d = NULL, max3d = NULL, showposition = TRUE, facets = NULL, plot = TRUE, storeData = TRUE, add = FALSE, pos = -1, ...) {
  ## dfxyz: matrix/data.frame with 3 columns
  ## min3d, max3d: limits by default: c(0,0,0), c(1,1,1)
  
  thecall <- .expand.call(match.call())
  
  ## parameters sorted
  sortparameters <- sortparamADEg(...)
  
  ## facets
  if(!is.null(facets)) {
    if(NCOL(fac) == 1)
      object <- multi.facets.Tr(thecall, samelimits = sortparameters$g.args$samelimits)
    else 
      stop("Facets are not allowed multiple fac")
  }
  
  ## multiple fac
  else if(NCOL(fac) > 1) {
    object <- multi.variables.Tr(thecall, "fac")
  }
  
  ## simple ADEg graphic
  else {
    if(length(sortparameters$rest))
      warning(c("Unused parameters: ", paste(unique(names(sortparameters$rest)), " ", sep = "")), call. = FALSE)
    
    ## creation of the ADEg object
    g.args <- c(sortparameters$g.args, list(adjust = adjust, min3d = min3d, max3d = max3d, col = col, order = thecall$order))
    if(storeData)
    	tmp_data <- list(dfxyz = dfxyz, fac = fac, labels = labels, frame = sys.nframe() + pos, storeData = storeData)
    else
      tmp_data <- list(dfxyz = thecall$dfxyz, fac = thecall$fac, labels = thecall$labels, frame = sys.nframe() + pos, storeData = storeData)
    object <- new(Class = "Tr.traject", data  = tmp_data, adeg.par = sortparameters$adepar, trellis.par = sortparameters$trellis, g.args = g.args, Call = match.call())
    
    ## preparation
    prepare(object)
    setlatticecall(object)
    if(showposition & add) {
      print("cannot show position and add") ## can be done, but modifies the meaning of the superposition
      showposition <- FALSE 
    }
    if(showposition)
      object <- new(Class = "ADEgS", ADEglist = list("triangle" = object, "positions" = .showpos(object)), positions = rbind(c(0, 0, 1, 1), c(0, 0.7, 0.3, 1)), add = matrix(0, ncol = 2, nrow = 2), Call = match.call())
    if(add)
      object <- add.ADEg(object)
  }
  
  if(!add & plot)
    print(object)
  invisible(object)
}

