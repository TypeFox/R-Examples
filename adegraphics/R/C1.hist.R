setClass(
  Class = "C1.hist",
  contains = "ADEg.C1"
)


setMethod(
  f = "initialize",
  signature = "C1.hist",
  definition = function(.Object, data = list(score = NULL, frame = 0, storeData = TRUE), ...) {
    .Object <- callNextMethod(.Object, data = data, ...) ## ADEg.C1 initialize
    validObject(.Object)
    return(.Object)
  })


setMethod(
  f = "prepare",
  signature = "C1.hist",
  definition = function(object) {
    nameobj <- deparse(substitute(object))
    
    ## pre-management of graphics parameters      
    oldparamadeg <- adegpar()
    on.exit(adegpar(oldparamadeg))
    adegtot <- adegpar(object@adeg.par)
    
    if(object@data$storeData)
      score <- object@data$score
    else
      score <- eval(object@data$score, envir = sys.frame(object@data$frame))
    
    score <- as.matrix(score)[, 1]  ## to manage 'score' when it is a data.frame with only one column
    
    ## change default for some parameters
    adegtot$p1d$rug$draw <- FALSE
    
    ## object modification before calling inherited method
    object@adeg.par <- adegtot
    callNextMethod() ## prepare graph
    
    ## compute histogram
    h <- hist(score, breaks = if(is.null(object@g.args$breaks)) object@g.args$nclass else object@g.args$breaks, right = object@g.args$right, plot = FALSE)
    y <- switch(object@g.args$type, count = h$counts, percent = 100 * h$counts / length(score), density = h$density)
    object@stats$heights <- y
    object@stats$breaks <- h$breaks
    
    if(object@adeg.par$p1d$horizontal && is.null(object@g.args$ylim))
	    object@g.args$ylim <- c(0, 1.1 * max(y))
    if(!object@adeg.par$p1d$horizontal && is.null(object@g.args$xlim))
	    object@g.args$xlim <- c(0, 1.1 * max(y))
    
    if(object@adeg.par$p1d$horizontal)
      object@g.args$scales$y$at <- pretty(object@g.args$ylim, n = 5)
    else
      object@g.args$scales$x$at <- pretty(object@g.args$xlim, n = 5)
    
    assign(nameobj, object, envir = parent.frame())
  })


setMethod(
  f = "panel",
  signature = "C1.hist",
  definition = function(object, x, y) {
    ## Drawing hist
    ## y is the score
    
    ## get some parameters    
    pscore <- object@adeg.par$p1d
    ppoly <- lapply(object@adeg.par$ppolygons, FUN = function(x) {rep(x, length.out = length(x))})
    breaks <- object@stats$breaks
    heights <- object@stats$heights
    
    ## Starts the display
    ## depends on the parametres horizontal
    ## reverse and rug.draw are always considered as FALSE
    if(pscore$horizontal) {
      panel.rect(x = breaks[-length(breaks)], y = 0, height = heights, width = diff(breaks), 
        col = ppoly$col, alpha = ppoly$alpha, border = ppoly$border, lty = ppoly$lty, 
        lwd = ppoly$lwd, just = c("left", "bottom"))
      
    } else {
      panel.rect(x = 0, y = breaks[-length(breaks)], height = diff(breaks), width = heights, 
        col = ppoly$col, alpha = ppoly$alpha, border = ppoly$border, lty = ppoly$lty, 
        lwd = ppoly$lwd, just = c("left", "bottom"))
    }
  })


s1d.hist <- function(score, breaks = NULL, nclass = round(log2(length(score)) + 1), type = c("count", "density", "percent"), right = TRUE, 
  facets = NULL, plot = TRUE, storeData = TRUE, add = FALSE, pos = -1, ...) {
  
  ## evaluation of some parameters
  thecall <- .expand.call(match.call())
  score <- eval(thecall$score, envir = sys.frame(sys.nframe() + pos))
  
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
    g.args <- c(sortparameters$g.args, list(type = match.arg(type), nclass = nclass, breaks = breaks, right = right))
    if(storeData)
    	tmp_data <- list(score = score, frame = sys.nframe() + pos, storeData = storeData)
    else
      tmp_data <- list(score = thecall$score, frame = sys.nframe() + pos, storeData = storeData)
    object <- new(Class = "C1.hist", data = tmp_data, adeg.par = sortparameters$adepar, trellis.par = sortparameters$trellis, g.args = g.args, Call = match.call())
    
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
