##########################################################################
##                            s.arrow                                   ##
##########################################################################

setClass(
  Class = "S2.arrow",
  contains = "ADEg.S2",
)


setMethod(
  f = "initialize",
  signature = "S2.arrow",
  definition = function(.Object, data = list(dfxy = NULL, xax = 1, yax = 2, labels = NULL, frame = 0, storeData = TRUE), ...) {
    .Object <- callNextMethod(.Object, data = data, ...) ## ADEg.S2 initialize
    .Object@data$labels <- data$labels
    return(.Object)
  })


setMethod(
  f = "prepare",
  signature = "S2.arrow",
  definition = function(object) {
    name_obj <- deparse(substitute(object))
    
    ## pre-management of graphics parameters      
    oldparamadeg <- adegpar()
    on.exit(adegpar(oldparamadeg))
    adegtot <- adegpar(object@adeg.par)
    
    ## change default for some parameters
    if(is.null(object@adeg.par$ppoints$cex))
      adegtot$ppoints$cex <- 0
    if(is.null(object@adeg.par$porigin$include) & (any(names(object@g.args) %in% c("Sp", "nbobject"))))
      adegtot$porigin$include <- FALSE
    
    ## object modification before calling inherited method
    object@adeg.par <- adegtot
    callNextMethod() ## prepare graph
    
    if(is.null(object@s.misc$lim.update)) {
      if(is.null(object@g.args$Sp)) {
        xdiff <- diff(object@g.args$xlim)
        ydiff <- diff(object@g.args$ylim)
        object@g.args$xlim <- object@g.args$xlim + c(-1, 1) * 0.05 * xdiff
        object@g.args$ylim <- object@g.args$ylim + c(-1, 1) * 0.05 * ydiff
      }
      object@s.misc$lim.update <- TRUE
    }
    assign(name_obj, object, envir = parent.frame())
  })


setMethod(
  f = "panel",
  signature = "S2.arrow",
  definition = function(object, x, y) {
    ## draw arrows
    panel.arrows(x0 = object@adeg.par$porigin$origin[1], y0 = object@adeg.par$porigin$origin[2], y1 = y, x1 = x, angle = object@adeg.par$parrows$angle, 
      					 length = object@adeg.par$parrows$length, ends = object@adeg.par$parrows$end, lwd = object@adeg.par$plines$lwd, 
                 col = object@adeg.par$plines$col, lty = object@adeg.par$plines$lty)
    ## draw labels
    ## positions
    plabels <- object@adeg.par$plabels
    if(object@data$storeData)
      arrownames <- object@data$labels
    else
      arrownames <- eval(object@data$labels, envir = sys.frame(object@data$frame))
    
    if(!is.null(arrownames)) {
      pos <- .textpos(x, y, origin = c(0, 0))
      test <- .textsize(arrownames, plabels)
      w <- test$w
      h <- test$h
      ## optim always false for s.arrow
      plabels$optim <- FALSE
      if(any(object@adeg.par$plabels$cex > 0))
        adeg.panel.label(x + pos[1, ] * w / 2, y + pos[2, ] * h / 2 , arrownames, plabels)
    }
  })


s.arrow <- function(dfxy, xax = 1, yax = 2, labels = row.names(as.data.frame(dfxy)), facets = NULL, plot = TRUE, storeData = TRUE, add = FALSE, pos = -1, ...) {
  
  ## evaluation of some parameters (required for multiplot)
  thecall <- .expand.call(match.call())
  df <- try(as.data.frame(eval(thecall$dfxy, envir = sys.frame(sys.nframe() + pos))), silent = TRUE)
  if((class(df) == "try-error") | is.null(thecall$dfxy)) ## non convenient dfxy argument
    stop("non convenient selection for dfxy (can not be converted to dataframe)")
  
  ## parameters sorted
  sortparameters <- sortparamADEg(...)

  ## facets
  if(!is.null(facets)) { 
    if((length(xax) == 1 & length(yax) == 1))
      object <- multi.facets.S2(thecall, sortparameters$adepar, samelimits = sortparameters$g.args$samelimits)
    else 
      stop("Facets are not allowed with multiple xax/yax")
  }
  
  ## multiple axes
  else if((length(xax) > 1 | length(yax) > 1)) {
    object <- multi.ax.S2(thecall)
  }
  
  ## simple ADEg graphic
  else {
    if(length(sortparameters$rest))
  	  warning(c("Unused parameters: ", paste(unique(names(sortparameters$rest)), " ", sep = "")), call. = FALSE)
    
    ## creation of the ADEg object
    if(storeData)
    	tmp_data <- list(dfxy = dfxy, xax = xax, yax = yax, labels = labels, frame = sys.nframe() + pos, storeData = storeData)
    else
      tmp_data <- list(dfxy = thecall$dfxy, xax = xax, yax = yax, labels = thecall$labels, frame = sys.nframe() + pos, storeData = storeData)
    object <- new(Class = "S2.arrow", data = tmp_data, adeg.par = sortparameters$adepar, trellis.par = sortparameters$trellis, g.args = sortparameters$g.args, Call = as.call(thecall))
    
    ## preparation
    prepare(object)
    setlatticecall(object)
    if(add)
      object <- add.ADEg(object)
  }

  if(!add & plot)
    print(object)
  invisible(object)  
}

