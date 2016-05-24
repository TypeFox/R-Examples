setClass(
  Class = "C1.dotplot",
  contains = "ADEg.C1"
)


setMethod(
  f = "initialize",
  signature = "C1.dotplot",
  definition = function(.Object, data = list(score = NULL, at = NULL, frame = 0, storeData = TRUE), ...) {
    .Object <- callNextMethod(.Object, data = data, ...) ## ADEg.C1 initialize
    .Object@data$at <- data$at
    validObject(.Object)
    return(.Object)
  })


setMethod(
  f = "prepare",
  signature = "C1.dotplot",
  definition = function(object) {
    nameobj <- deparse(substitute(object))
    
    ## pre-management of graphics parameters      
    oldparamadeg <- adegpar()
    on.exit(adegpar(oldparamadeg))
    adegtot <- adegpar(object@adeg.par)
    
    if(object@data$storeData)
      at <- object@data$at
    else
      at <- eval(object@data$at, envir = sys.frame(object@data$frame))
    
    ## change some defaults
    adegtot$p1d$rug$draw <- FALSE
    ## object modification before calling inherited method
    object@adeg.par <- adegtot
    callNextMethod() ## prepare graph
    
    if(object@adeg.par$p1d$horizontal && is.null(object@g.args$ylim))
      object@g.args$ylim <- setlimits1D(min(at), max(at), 0, FALSE)  
    if(!object@adeg.par$p1d$horizontal && is.null(object@g.args$xlim))
      object@g.args$xlim <- setlimits1D(min(at), max(at), 0, FALSE)
    
    assign(nameobj, object, envir = parent.frame())
  })


setMethod(
  f = "panel",
  signature = "C1.dotplot",
  definition = function(object, x, y) {
    ## Drawing dotchart
    ## x is the index
    ## y is the score
    
    ## get some parameters    
    pscore <- object@adeg.par$p1d
    ppoints <- lapply(object@adeg.par$ppoints, FUN = function(x) {rep(x, length.out = length(x))})
    plines <- lapply(object@adeg.par$plines, FUN = function(x) {rep(x, length.out = length(x))})
    
    ## reorder the values
    y <- y[order(x)]
    x <- sort(x)
    
    ## Starts the display
    ## depends on the parametres horizontal
    ## rug.draw and reverse are always considered as FALSE
    
    if(pscore$horizontal) {
      x.tmp <- y
      y.tmp <- x
      panel.segments(object@adeg.par$porigin$origin[1], y.tmp, x.tmp, y.tmp, lwd = plines$lwd, lty = plines$lty, col = plines$col)
    } else {
      x.tmp <- x
      y.tmp <- y
      panel.segments(x.tmp, object@adeg.par$porigin$origin[1], x.tmp, y.tmp, lwd = plines$lwd, lty = plines$lty, col = plines$col)
    }
    
    panel.dotplot(x = x.tmp, y = y.tmp, horizontal = pscore$horizontal, pch = ppoints$pch, cex = ppoints$cex, col = ppoints$col, alpha = ppoints$alpha, col.line = "transparent")
  })


s1d.dotplot <- function(score, at = 1:NROW(score), facets = NULL, plot = TRUE, storeData = TRUE, add = FALSE, pos = -1, ...) {
  
  ## evaluation of some parameters
  thecall <- .expand.call(match.call())
  ## parameters sorted
  sortparameters <- sortparamADEg(...)
  
  ## facets
  if(!is.null(facets)) {
    if(NCOL(score) == 1)
      object <- multi.facets.C1(thecall, sortparameters$adepar, samelimits = sortparameters$g.args$samelimits)
    else 
      stop("Facets are not allowed with multiple scores")
  }
  
  ## multiple scores
  else if(NCOL(score) > 1) {
    object <- multi.score.C1(thecall)
  }
  
  ## simple ADEg graphic
  else {
    if(length(sortparameters$rest))
      warning(c("Unused parameters: ", paste(unique(names(sortparameters$rest)), " ", sep = "")), call. = FALSE)
    
    ## creation of the ADEg object
    if(storeData)
    	tmp_data <- list(score = score, at = at, frame = sys.nframe() + pos, storeData = storeData)
    else
      tmp_data <- list(score = thecall$score, at = thecall$at, frame = sys.nframe() + pos, storeData = storeData)
    object <- new(Class = "C1.dotplot", data = tmp_data, adeg.par = sortparameters$adepar, trellis.par = sortparameters$trellis, g.args = sortparameters$g.args, Call = match.call())
    
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
