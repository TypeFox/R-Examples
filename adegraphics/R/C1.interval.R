setClass(
  Class = "C1.interval",
  contains = "ADEg.C1"
)


setMethod(
  f = "initialize",
  signature = "C1.interval",
  definition = function(.Object, data = list(score = NULL, at = NULL, frame = 0, storeData = TRUE), ...) {
    .Object <- callNextMethod(.Object, data = data, ...) ## ADEg.C1 initialize
    .Object@data$at <- data$at
    validObject(.Object)
    return(.Object)
  })


setMethod(
  f = "prepare",
  signature = "C1.interval",
  definition = function(object) {
    name_obj <- deparse(substitute(object))
    
    ## pre-management of graphics parameters      
    oldparamadeg <- adegpar()
    on.exit(adegpar(oldparamadeg))
    adegtot <- adegpar(object@adeg.par)
    
    if(object@data$storeData)
      at <- object@data$at
    else
      at <- eval(object@data$at, envir = sys.frame(object@data$frame))
    
    ## change default for some parameters
    adegtot$p1d$rug$draw <- FALSE
    if(object@g.args$method == "bars") {
      if(is.null(object@adeg.par$parrows$ends))
        adegtot$parrows$ends <- "both"
      if(is.null(object@adeg.par$parrows$angle))
        adegtot$parrows$angle <- 90    
    }     
    
    ## object modification before calling inherited method
    object@adeg.par <- adegtot
    callNextMethod() ## prepare graph
    
    if(object@adeg.par$p1d$horizontal && is.null(object@g.args$ylim))
      object@g.args$ylim <- setlimits1D(min(at), max(at), 0, FALSE)
    if(!object@adeg.par$p1d$horizontal && is.null(object@g.args$xlim))
      object@g.args$xlim <- setlimits1D(min(at), max(at), 0, FALSE)
    
    assign(name_obj, object, envir = parent.frame())
  })


setMethod(
  f= "panel",
  signature = "C1.interval",
  definition = function(object, x, y) {
    ## Drawing interval
    ## x is the index
    ## y is the score  
    
    lims <- current.panel.limits(unit = "native")
    
    pscore <- object@adeg.par$p1d
    plines <- object@adeg.par$plines
    parrows <- object@adeg.par$parrows
    ppoly <- object@adeg.par$ppolygons
    
    nval <- length(y) %/% 2
    score2 <- y[(nval + 1):length(y)]
    score1 <- y[1 : nval]
    
    ## reorder the values
    score1 <- score1[order(x)]
    score2 <- score2[order(x)]
    x <- sort(x)
    
    ## Starts the display
    ## depends on the parametres horizontal
    ## rug.draw and reverse are always considered as FALSE
    
    if(pscore$horizontal) {
      if(object@g.args$method == "area") {
        panel.polygon(x = c(score1, rev(score2)), y = c(x, rev(x)), border = "transparent", col = ppoly$col, alpha = ppoly$alpha)
        panel.lines(x = score1, y = x, col = ppoly$border, lty = ppoly$lty, lwd = ppoly$lwd)
        panel.lines(x = score2, y = x, col = ppoly$border, lty = ppoly$lty, lwd = ppoly$lwd)
      } else if(object@g.args$method == "bars") {
        panel.arrows(x0 = score1, y0 = x, x1 = score2, y1 = x, lwd = plines$lwd, col = plines$col, 
          lty = plines$lty, angle = parrows$angle, length = parrows$length, ends = parrows$ends)
      }
      
    } else {
      if(object@g.args$method == "area") {
        panel.polygon(x = c(x, rev(x)), y = c(score1, rev(score2)),  border = "transparent", col = ppoly$col, alpha = ppoly$alpha)
        panel.lines(x = x, y = score1, col = ppoly$border, lty = ppoly$lty, lwd = ppoly$lwd)
        panel.lines(x = x, y = score2, col = ppoly$border, lty = ppoly$lty, lwd = ppoly$lwd)
      } else if(object@g.args$method == "bars") {
        panel.arrows(x0 = x, y0 = score1, x1 = x, y1 = score2, lwd = plines$lwd, col = plines$col, 
          lty = plines$lty, angle = parrows$angle, length = parrows$length, ends = parrows$ends)
      }
    }
  })


s1d.interval <- function(score1, score2, at = 1:NROW(score1), method = c("bars", "area"), facets = NULL, plot = TRUE, storeData = TRUE, add = FALSE, pos = -1, ...) {
  
  ## evaluation of some parameters
  thecall <- .expand.call(match.call())
  score1 <- eval(thecall$score1, envir = sys.frame(sys.nframe() + pos))
  score2 <- eval(thecall$score2, envir = sys.frame(sys.nframe() + pos))
  if(NROW(score1) != NROW(score2))
    stop("score1 and score2 should have the same length")
  if(NCOL(score1) != NCOL(score2))
    stop("score1 and score2 should have the same number of columns")
  
  if((is.data.frame(score1) & NCOL(score1) == 1) | (is.data.frame(score2) & NCOL(score2) == 1)) 
    stop("Not yet implemented for data.frame with only one column, please convert into vector")
  
  ## parameters sorted
  sortparameters <- sortparamADEg(...)
  
  ## facets
  if(!is.null(facets)) {
    if(NCOL(score1) == 1)
      object <- multi.facets.C1(thecall, sortparameters$adepar, samelimits = sortparameters$g.args$samelimits)
    else 
      stop("Facets are not allowed with multiple scores")
  }
  
  ## multiple scores
  else if(NCOL(score1) > 1) { 
    object <- multi.score.C1(thecall)
  }
  
  ## simple ADEg graphic
  else {
    if(length(sortparameters$rest))
      warning(c("Unused parameters: ", paste(unique(names(sortparameters$rest)), " ", sep = "")), call. = FALSE)
    
    ## creation of the ADEg object
    g.args <- c(sortparameters$g.args, list(method = match.arg(method)))
    if(storeData)
    	tmp_data <- list(score = c(score1, score2), at = at, frame = sys.nframe() + pos, storeData = storeData)
    else
      tmp_data <- list(score = call("c", thecall$score1, thecall$score2), at = thecall$at, frame = sys.nframe() + pos, storeData = storeData)
    object <- new(Class = "C1.interval", data = tmp_data, adeg.par = sortparameters$adepar, trellis.par = sortparameters$trellis, g.args = g.args, Call = match.call())
    
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
