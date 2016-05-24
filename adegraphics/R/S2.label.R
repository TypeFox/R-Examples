#########################################################
###                       s.label                      ##
#########################################################

setClass(
  Class = "S2.label",
  contains = "ADEg.S2"
)


setMethod(
  f = "initialize",
  signature = "S2.label",
  definition = function(.Object, data = list(dfxy = NULL, labels = NULL, xax = 1, yax = 2, frame = 0, storeData = TRUE), ...) {
    .Object <- callNextMethod(.Object, data = data, ...) ## ADEg.S2 initialize
    .Object@data$labels <- data$labels
    return(.Object)
  })


setMethod(
  f = "prepare",
  signature = "S2.label",
  definition = function(object) {
    name_obj <- deparse(substitute(object))
    
    if(object@data$storeData) {
      labels <- object@data$labels
    } else {
      labels <- eval(object@data$labels, envir = sys.frame(object@data$frame))
    }
    
    ## pre-management of graphics parameters      
    oldparamadeg <- adegpar()
    on.exit(adegpar(oldparamadeg))
    adegtot <- adegpar(object@adeg.par)
    
    ## change default for some parameters
    if((is.null(object@adeg.par$plabels$boxes$draw) & adegtot$plabels$optim) || (is.null(object@adeg.par$plabels$boxes$draw) & length(labels) > 1000)) ## dessin des boxes labels non renseigne par l'utilisateur lors de lappel a factorise
      adegtot$plabels$boxes$draw <- FALSE
    if(is.null(object@adeg.par$porigin$include) & (any(names(object@g.args) %in% c("Sp", "nbobject"))))
      adegtot$porigin$include <- FALSE
    
    ## object modification before calling inherited method
    object@adeg.par <- adegtot
    callNextMethod() ## prepare graph
    assign(name_obj, object, envir = parent.frame())
  })


setMethod(
  f = "panel",
  signature = "S2.label",
  definition = function(object, x, y) {
    ## draw labels
    if(any(object@adeg.par$ppoints$cex > 0))
      panel.points(x, y, pch = object@adeg.par$ppoints$pch, cex = object@adeg.par$ppoints$cex, col = object@adeg.par$ppoints$col, alpha = object@adeg.par$ppoints$alpha, fill = object@adeg.par$ppoints$fill)
    
    if(object@data$storeData)
      labels <- object@data$labels
    else
      labels <- eval(object@data$labels, envir = sys.frame(object@data$frame))
    
    if(any(object@adeg.par$plabels$cex > 0) & (!is.null(labels)))
      adeg.panel.label(x, y, labels, object@adeg.par$plabels)
  })


s.label <- function(dfxy, labels = rownames(dfxy), xax = 1, yax = 2, facets = NULL, plot = TRUE, storeData = TRUE, add = FALSE, pos = -1, ...) {
  
  ## evaluation of some parameters
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
    object <- new(Class = "S2.label", data = tmp_data, adeg.par = sortparameters$adepar, trellis.par = sortparameters$trellis, g.args = sortparameters$g.args, Call = as.call(thecall))
    
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
