#########################################################
###                       s.match                      ##
#########################################################

## in S2.match, the two data set are combined (using rbind) and kept as the same one...
## We know that the two data sets have the same row number, so we can easily retrieve and distinguish the two set (the first (nrow/2) rows are from dfxy1 the rest from dfxy2
setClass(
  Class = "S2.match",
  contains = "ADEg.S2",
)


setMethod(
  f = "initialize",
  signature = "S2.match",
  definition = function(.Object, data = list(dfxy = NULL, xax = 1, yax = 2, labels = NULL, frame = 0, storeData = TRUE), ...) {
    .Object <- callNextMethod(.Object, data = data, ...) ## ADEg.S2 initialize
    .Object@data$labels <- data$labels
    return(.Object)
  })


setMethod(
  f = "prepare",
  signature = "S2.match",
  definition = function(object) {
    name_obj <- deparse(substitute(object))
    
    ## pre-management of graphics parameters      
    oldparamadeg <- adegpar()
    on.exit(adegpar(oldparamadeg))
    adegtot <- adegpar(object@adeg.par)
    
    ## change default for some parameters
    if(!object@g.args$arrows)
      adegtot$parrows$angle <- 0
    else
    	if(is.null(object@adeg.par$parrows$angle) || (object@adeg.par$parrows$angle == 0))
	      adegtot$parrows$angle <- oldparamadeg$parrows$angle
      
    if(is.null(object@adeg.par$porigin$include) & (any(names(object@g.args) %in% c("Sp", "nbobject"))))
      adegtot$porigin$include <- FALSE
    
    ## object modification before calling inherited method
    object@adeg.par <- adegtot
    callNextMethod() ## prepare graph
    assign(name_obj, object, envir = parent.frame())
  })


setMethod(
  f = "panel",
  signature = "S2.match",
  definition = function(object, x, y) {
    n <- length(x) / 2
    if(length(x) %% 2) ## x non multiple de 2
      stop("error in spanel, not finding two datasets with equal row numbers please see with dev")
    ## arrows from dfxy to dfxy2
    panel.arrows(x0 = x[1:n], y0 = y[1:n] , y1 = y[n + 1:(2 * n)], x1 = x[n + 1:(2 * n)], angle = object@adeg.par$parrows$angle,
                 length = object@adeg.par$parrows$length, ends = object@adeg.par$parrows$end,
                 lwd = object@adeg.par$plines$lwd, col = object@adeg.par$plines$col, lty = object@adeg.par$plines$lty)

    do.call("panel.points", c(list(x = x[1:n], y = y[1:n]), object@adeg.par$ppoints))
    ## dessins labels
    if(any(object@adeg.par$plabels$cex > 0)) {
      xlab <- ((x[1:n] + x[(n + 1):(2 * n)]) / 2)
      ylab <- ((y[1:n] + y[(n + 1):(2 * n)]) / 2)
      if(object@data$storeData)
        labels <- object@data$labels
      else
        labels <- eval(object@data$labels, envir = sys.frame(object@data$frame))
      object@adeg.par$plabels$optim <- FALSE
      adeg.panel.label(xlab, ylab ,labels, object@adeg.par$plabels)
    }
  })


## if arrows= TRUE arrows are plotted, otherwise only the segments are drawn
s.match <- function(dfxy1, dfxy2, xax = 1, yax = 2, labels = row.names(as.data.frame(dfxy1)), arrows = TRUE, facets = NULL, plot = TRUE,
  storeData = TRUE, add = FALSE, pos = -1, ...) {

  ## evaluation of some parameters
  thecall <- .expand.call(match.call())
  data1 <- try(as.data.frame(eval(thecall$dfxy1, envir = sys.frame(sys.nframe() + pos))), silent = TRUE)
  data2 <- try(as.data.frame(eval(thecall$dfxy2, envir = sys.frame(sys.nframe() + pos))), silent = TRUE)
  
  if(class(data1) == "try-error" || class(data2) == "try-error" || is.null(thecall$dfxy1) || is.null(thecall$dfxy2))  ## wrong conversion 
    stop("non convenient selection for dfxy1 or dfxy2 (can not be converted to dataframe)")
  if(any(is.na(pmatch(colnames(data1), colnames(data2)))))
    stop("column names should be identical")
  if(any(is.na(data1)))
    stop("NA in first dataframe") ## TODO
  if(any(is.na(data2)))
    stop("NA in second dataframe") ## TODO
  if(nrow(data1) != nrow(data2)) 
    stop("non equal row numbers")
  
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
    g.args <- c(sortparameters$g.args, list(arrows = arrows))
    if(storeData)
    	tmp_data <- list(dfxy = rbind(dfxy1, dfxy2), xax = xax, yax = yax, labels = labels, frame = sys.nframe() + pos, storeData = storeData)
    else
      tmp_data <- list(dfxy = call("rbind", thecall$dfxy1, thecall$dfxy2), xax = xax, yax = yax, labels = thecall$labels, frame = sys.nframe() + pos, storeData = storeData)
    object <- new(Class = "S2.match", data = tmp_data , adeg.par = sortparameters$adepar, trellis.par = sortparameters$trellis, g.args = g.args, Call = as.call(thecall))
    
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

