#########################################################
###                     s.traject                      ##
#########################################################

setClass(
  Class= "S2.traject",
  contains = "ADEg.S2",
)


setMethod(
  f = "initialize",
  signature = "S2.traject",
  definition = function(.Object, data = list(dfxy = NULL, fac = NULL, labels = NULL, xax = 1, yax = 2, frame = 0, storeData = TRUE), ...) {
    .Object <- callNextMethod(.Object, data = data, ...) ## ADEg.S2 initialize
    .Object@data$labels <- data$labels
    .Object@data$fac <- data$fac
    return(.Object)
  })


setMethod(
  f = "prepare",
  signature = "S2.traject",
  definition = function(object) {
    name_obj <- deparse(substitute(object))
    
    if(object@data$storeData)
      fac <- as.factor(object@data$fac)
    else
      fac <- as.factor(eval(object@data$fac, envir = sys.frame(object@data$frame)))
    
    ## pre-management of graphics parameters      
    oldparamadeg <- adegpar()
    on.exit(adegpar(oldparamadeg))
    adegtot <- adegpar(object@adeg.par)
    
    ## matching col, pch ad size for points
    if(!is.null(object@g.args$col))
      if(is.logical(object@g.args$col)) {
        if(object@g.args$col)
          adegtot$ppoints$col <- adegtot$ppoints$fill <- adegtot$plabels$col <- adegtot$plabels$boxes$border <- adegtot$plines$col <- adegtot$ppalette$quali(nlevels(fac))
      } else
        adegtot$ppoints$col <- adegtot$ppoints$fill <- adegtot$plabels$col <- adegtot$plabels$boxes$border <- adegtot$plines$col <- rep(object@g.args$col, length.out = nlevels(fac))
    
    if(is.null(object@adeg.par$porigin$include) & (any(names(object@g.args) %in% c("Sp", "nbobject"))))
      adegtot$porigin$include <- FALSE
    
    ## object modification before calling inherited method
    object@adeg.par <- adegtot
    callNextMethod() ## prepare graph
    assign(name_obj, object, envir = parent.frame())
  })


setMethod(
  f = "panel",
  signature = "S2.traject",
  definition = function(object, x, y) {
    if(object@data$storeData) {
      fact <- object@data$fac
      labels <- object@data$labels
    } else {
      fact <- eval(object@data$fac, envir = sys.frame(object@data$frame))
      labels <- eval(object@data$labels, envir = sys.frame(object@data$frame))
    }
    
    todrawX <- split(x, fact)
    todrawY <- split(y, fact)
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


s.traject <- function(dfxy, fac = gl(1, nrow(dfxy)), order, labels = levels(fac), xax = 1, yax = 2, col = NULL, facets = NULL, plot = TRUE, storeData = TRUE, add = FALSE, pos = -1, ...) {
  
  ## evaluation of some parameters (required for multiplot)
  thecall <- .expand.call(match.call())
  labels <- eval(thecall$labels, envir = sys.frame(sys.nframe() + pos))
  fac <- eval(thecall$fac, envir = sys.frame(sys.nframe() + pos))
  df <- try(as.data.frame(eval(thecall$dfxy, envir = sys.frame(sys.nframe() + pos))), silent = TRUE)
  if((class(df) == "try-error") | is.null(thecall$dfxy)) ## non convenient dfxy argument
    stop("non convenient selection for dfxy (can not be converted to dataframe)")
  
  if(missing(fac))
    stop("no factor specified")
  
  if(NCOL(fac) == 1) {
    fac <- as.factor(fac)
    if(length(labels) != nlevels(fac))
      stop("wrong number of labels")
  }
  
  ## parameters sorted
  sortparameters <- sortparamADEg(...)
  
  ## facets
  if(!is.null(facets)) { 
    if((length(xax) == 1 & length(yax) == 1) & NCOL(fac) == 1)
      object <- multi.facets.S2(thecall, sortparameters$adepar, samelimits = sortparameters$g.args$samelimits)
    else 
      stop("Facets are not allowed with multiple xax/yax or multiple fac")
  }
  
  ## multiple axes
  else if((length(xax) > 1 | length(yax) > 1)) {
    if(NCOL(fac) == 1)
      object <- multi.ax.S2(thecall)
    else 
      stop("Multiple xax/yax are not allowed with multiple fac")
  }
  
  ## multiple fac
  else if(NCOL(fac) > 1) {
    object <- multi.variables.S2(thecall, "fac")
  }
  
  ## simple ADEg graphic
  else {
    if(length(sortparameters$rest))
      warning(c("Unused parameters: ", paste(unique(names(sortparameters$rest)), " ", sep = "")), call. = FALSE)
    
    ## creation of the ADEg object
    g.args <- c(sortparameters$g.args, list(order = thecall$order, col = col))
    if(storeData)
    	tmp_data <- list(dfxy = dfxy, xax = xax, yax = yax, labels = labels, fac = fac, frame = sys.nframe() + pos, storeData = storeData)
    else
      tmp_data <- list(dfxy = thecall$dfxy, xax = xax, yax = yax, labels = thecall$labels, fac = thecall$fac, frame = sys.nframe() + pos, storeData = storeData)
    object <- new(Class = "S2.traject", data = tmp_data, adeg.par = sortparameters$adepar, trellis.par = sortparameters$trellis, g.args = g.args, Call = as.call(thecall))
    
    ## preparation
    prepare(object)
    setlatticecall(object)
    if(add)
      object <- add.ADEg(object)
  }
  
  if(! add & plot)
    print(object)
  invisible(object)
}
