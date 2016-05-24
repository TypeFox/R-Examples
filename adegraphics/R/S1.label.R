###########################################################
##                          s1d.label                    ##
###########################################################

setClass(
  Class = "S1.label",
  contains = "ADEg.S1"
)


setMethod(
  f = "initialize",
  signature = "S1.label",
  definition = function(.Object, data = list(score = NULL, labels = NULL, at = NULL, frame = 0, storeData = TRUE), ...) {
    .Object <- callNextMethod(.Object, data = data, ...) ## ADEg.C1 initialize
    .Object@data$labels <- data$labels
    return(.Object)
  })


setMethod(
  f = "prepare",
  signature = "S1.label",
  definition = function(object) {
    name_obj <- deparse(substitute(object))
    
    ## pre-management of graphics parameters      
    oldparamadeg <- adegpar()
    on.exit(adegpar(oldparamadeg))
    adegtot <- adegpar(object@adeg.par)
    
    ## change default for some parameters
    if(adegtot$p1d$horizontal & is.null(object@adeg.par$plabels$srt))
      adegtot$plabels$srt <- 90
    else if(!adegtot$p1d$horizontal & is.null(object@adeg.par$plabels$srt))
      adegtot$plabels$srt <- 0
    
    if(adegtot$p1d$horizontal & is.null(object@g.args$ylim))
      object@g.args$ylim <- c(0, 1)
    
    if(!adegtot$p1d$horizontal & is.null(object@g.args$xlim))
      object@g.args$xlim <- c(0, 1)
    
    ## object modification before calling inherited method
    object@adeg.par <- adegtot
    callNextMethod() ## prepare graph
    
    assign(name_obj, object, envir = parent.frame())
  })


setMethod(
  f = "panel",
  signature = "S1.label",
  definition = function(object, x, y) {
    
    if(object@data$storeData) {
      labels <- object@data$labels
      at <- object@data$at
    } else {
      labels <- eval(object@data$labels, envir = sys.frame(object@data$frame))
      at <- eval(object@data$at, envir = sys.frame(object@data$frame))
    }
    
    lims <- current.panel.limits(unit = "native")
    pscore <- object@adeg.par$p1d
    plabels <- object@adeg.par$plabels
    plboxes <- plabels$boxes
    nval <- length(y)
    
    if(!is.null(labels)) {
      ## get text sizes for boxes
      test <- .textsize(labels, plabels)
      w <- test$w
      h <- test$h
    }        
    
    lead <- ifelse(pscore$reverse, -1, 1)
    
    if(pscore$horizontal) {
      ## horizontal plot
      xpoints <- y
      
      ## draw labels
      if(object@g.args$poslabel == "regular") {
        spacelab <- diff(lims$xlim) / (nval + 1)
        xlab <- seq(from = lims$xlim[1] + spacelab, by = spacelab, length.out = nval)[rank(xpoints, ties.method = "first")]
      } else
        xlab <- xpoints
      
      if(!is.null(labels) & any(plabels$cex > 0))
        adeg.panel.label(x = xlab , y = at + lead * h / 2, labels = labels, plabels = plabels)
      
      ## draw segments
      ypoints <- object@s.misc$rug
      do.call("panel.segments", c(list(x0 = xpoints, y0 = ypoints, x1 = xlab, y1 = at), object@adeg.par$plines))
      
    } else {
      ## vertical plot
      ypoints <- y
      
      ## draw labels
      if(object@g.args$poslabel == "regular") {
        spacelab <- diff(lims$ylim) / (nval + 1)
        ylab <- seq(from = lims$ylim[1] + spacelab, by = spacelab, length.out = nval)[rank(ypoints, ties.method = "first")]
      } else
        ylab <- ypoints
      if(!is.null(labels) & any(plabels$cex > 0))
        adeg.panel.label(x = at + lead * w / 2 , y = ylab, labels = labels, plabels = plabels)
      
      ## draw segments
      xpoints <- object@s.misc$rug
      do.call("panel.segments", c(list(x0 = xpoints, y0 = ypoints, x1 = at, y1 = ylab), object@adeg.par$plines))
    }
    
    if(any(object@adeg.par$ppoints$cex > 0))
      do.call("panel.points", c(list(x = xpoints, y = ypoints), object@adeg.par$ppoints))
  })


s1d.label <- function(score, labels = 1:NROW(score), at = 0.5, poslabel = c("regular", "value"), facets = NULL, plot = TRUE, storeData = TRUE, add = FALSE, pos = -1, ...) {
  
  ## evaluation of some parameters
  thecall <- .expand.call(match.call())
  score <- eval(thecall$score, envir = sys.frame(sys.nframe() + pos))
  
  ## parameters sorted
  sortparameters <- sortparamADEg(...)
  
  ## facets
  if(!is.null(facets)) {
    if(NCOL(score) == 1)
      object <- multi.facets.S1(thecall, sortparameters$adepar, samelimits = sortparameters$g.args$samelimits)
    else 
      stop("Facets are not allowed with multiple scores")
  }
  
  ## multiple scores
  else if(NCOL(score) > 1) { 
    object <- multi.score.S1(thecall)
  }
  
  ## simple ADEg graphic
  else {
    if(length(sortparameters$rest))
      warning(c("Unused parameters: ", paste(unique(names(sortparameters$rest)), " ", sep = "")), call. = FALSE)
    
    ## creation of the ADEg object
    g.args <- c(sortparameters$g.args, list(poslabel = match.arg(poslabel)))
    if(storeData)
      tmp_data <- list(score = score, labels = labels, at = at, frame = sys.nframe() + pos, storeData = storeData)
    else
      tmp_data <- list(score = thecall$score, labels = thecall$labels, at = thecall$at, frame = sys.nframe() + pos, storeData = storeData)
    object <- new(Class = "S1.label", data = tmp_data, adeg.par = sortparameters$adepar, trellis.par = sortparameters$trellis, g.args = g.args, Call = match.call())
    
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
