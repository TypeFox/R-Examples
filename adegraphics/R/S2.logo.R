#########################################################
##                         s.logo                      ##
#########################################################


setClass(
  Class = "S2.logo",
  contains = "ADEg.S2",
)

setMethod(
  f = "initialize",
  signature = "S2.logo",
  definition = function(.Object, data = list(dfxy = NULL, logos = NULL, xax = 1, yax = 2, frame = 0, storeData = TRUE), ...) {
    .Object <- callNextMethod(.Object, data = data, ...) ## ADEg.S2 initialize
    .Object@data$logos <- data$logos
    return(.Object)
  })


setMethod(
  f = "prepare",
  signature = "S2.logo",
  definition = function(object) {
    name_obj <- deparse(substitute(object))
    
    ## pre-management of graphics parameters      
    oldparamadeg <- adegpar()
    on.exit(adegpar(oldparamadeg))
    adegtot <- adegpar(object@adeg.par)
    
    if(is.null(object@adeg.par$porigin$include) & (any(names(object@g.args) %in% c("Sp", "nbobject"))))
      adegtot$porigin$include <- FALSE
    
    ## object modification before calling inherited method
    object@adeg.par <- adegtot
    callNextMethod() ## prepare graph
    assign(name_obj, object, envir = parent.frame())
  })


setMethod(
  f = "panel",
  signature = "S2.logo",
  definition = function(object, x, y) {
    ## list of bitmap objects:logos
    if(object@data$storeData)
      logos <- object@data$logos
    else
      logos <- eval(object@data$logos, envir = sys.frame(object@data$frame))
    
    for(i in 1:length(logos)) {
      grid.draw(rasterGrob(logos[[i]], x = x[i], y = y[i], height = unit(0.1, "npc") * object@adeg.par$ppoints$cex, default.units = "native"))
    }
  })


s.logo <- function(dfxy, logos, xax = 1, yax = 2, facets = NULL, plot = TRUE, storeData = TRUE, add = FALSE, pos = -1, ...) {
  
  ## evaluation of some parameters
  thecall <- .expand.call(match.call())
  df <- try(as.data.frame(eval(thecall$dfxy, envir = sys.frame(sys.nframe() + pos))), silent = TRUE)
  if((class(df) == "try-error") | is.null(thecall$dfxy)) ## non convenient dfxy argument
    stop("non convenient selection for dfxy (can not be converted to dataframe)")
  logos <- eval(thecall$logos, envir = sys.frame(sys.nframe() + pos))
  if(!is.list(logos))
    stop("The argument 'logos' should be a list")
  
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
    	tmp_data <- list(dfxy = dfxy, xax = xax, yax = yax, logos = logos, frame = sys.nframe() + pos, storeData = storeData)
    else
      tmp_data <- list(dfxy = thecall$dfxy, xax = xax, yax = yax, logos = thecall$logos, frame = sys.nframe() + pos, storeData = storeData)
    object <- new(Class = "S2.logo", data = tmp_data, adeg.par = sortparameters$adepar, trellis.par = sortparameters$trellis, g.args = sortparameters$g.args, Call =  as.call(thecall))
    
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

