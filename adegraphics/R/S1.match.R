###########################################################
##                          s1d.match                    ##
###########################################################

setClass(
  Class = "S1.match",
  contains = "ADEg.S1"
)


setMethod(
  f = "initialize",
  signature = "S1.match",
  definition = function(.Object, data = list(score = NULL, labels = NULL, at = NULL, frame = 0, storeData = TRUE), ...) {
    .Object <- callNextMethod(.Object, data = data, ...) ## ADEg.S1 initialize
    .Object@data$labels <- data$labels
    return(.Object)
  })


setMethod(
  f = "prepare",
  signature = "S1.match",
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
    adegtot$p1d$rug$tck <- 0
    
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
  f= "panel",
  signature = "S1.match",
  definition = function(object, x, y) {
    
    if(object@data$storeData) {
      labels <- object@data$labels
      at <- object@data$at
    } else {
      labels <- eval(object@data$labels, envir = sys.frame(object@data$frame))
      at <- eval(object@data$at, envir = sys.frame(object@data$frame))
    }
    
    lims <- current.panel.limits(unit = "native")
    nval <- length(y) %/% 2
    score2 <- y[(nval + 1):length(y)]
    score1 <- y[1 : nval]
    
    pscore <- object@adeg.par$p1d
    plabels <- object@adeg.par$plabels
    plboxes <- plabels$boxes
    porigin <- object@adeg.par$porigin
    
    if(!is.null(labels)) {
      ## get text sizes for boxes
      test <- .textsize(labels, plabels)
      w <- test$w
      h <- test$h
    }
    
    lead <- ifelse(pscore$reverse, -1, 1)
    
    if(pscore$horizontal) {
      ## horizontal plot
      ## get positions for labels
      spacelab <- diff(lims$xlim) / (nval + 1)
      xlab <- seq(from = lims$xlim[1] + spacelab, by = spacelab, length.out = nval)[rank(score1, ties.method = "first")]
      ylab <- rep(at, length.out = nval)
      
      ypoints <- rep(object@s.misc$rug, length.out = nval)
      ypoints2 <- rep(ypoints + lead * 0.05 * abs(diff(object@g.args$ylim)), length.out = nval)
      
      ## horizontal line
      if(pscore$rug$draw & pscore$rug$line) 
        panel.abline(h = ypoints2, col = porigin$col, lwd = porigin$lwd, lty = porigin$lty, alpha = porigin$alpha)
      ## segments linking both scores
      do.call("panel.segments", c(list(x0 = score1, y0 = ypoints, x1 = score2, y1 = ypoints2), object@adeg.par$plines))
      ## segments linking labels to second score
      do.call("panel.segments", c(list(x0 = score2, y0 = ypoints2, x1 = xlab, y1 = ylab), object@adeg.par$plines))
      
      ## drawing labels
      if(!is.null(labels) & any(plabels$cex > 0))
        adeg.panel.label(x = xlab , y = ylab + lead * h / 2, labels = labels, plabels = plabels)
      ## draw points
      if(any(object@adeg.par$ppoints$cex > 0))
        panel.points(x = c(score1, score2), y = c(ypoints, ypoints2), pch = object@adeg.par$ppoints$pch, cex = object@adeg.par$ppoints$cex, col = object@adeg.par$ppoints$col, alpha = object@adeg.par$ppoints$alpha, fill = object@adeg.par$ppoints$fill)
      
    } else {
      ## vertical plot
      ## get positions for labels
      spacelab <- diff(lims$ylim) / (nval + 1)
      ylab <- seq(from = lims$ylim[1] + spacelab, by = spacelab, length.out = nval)[rank(score1, ties.method = "first")]
      xlab <- rep(at, length.out = nval)
      
      xpoints <- rep(object@s.misc$rug, length.out = nval)
      xpoints2 <- rep(xpoints + lead * 0.05 * abs(diff(object@g.args$xlim)), length.out = nval)
      
      ## vertical line
      if(pscore$rug$draw & pscore$rug$line) 
        panel.abline(v = xpoints2,  col = porigin$col, lwd = porigin$lwd, lty = porigin$lty, alpha = porigin$alpha)
      ## segments linking both scores
      do.call("panel.segments", c(list(x0 = xpoints, y0 = score1, x1 =  xpoints2, y1 = score2), object@adeg.par$plines))
      ## segments linking labels to second score
      do.call("panel.segments", c(list(x0 = xpoints2, y0 = score2, x1 = xlab, y1 = ylab), object@adeg.par$plines))
      
      ## drawing labels
      if(!is.null(labels) & any(plabels$cex > 0))
        adeg.panel.label(x = xlab + lead * w / 2 , y = ylab, labels = labels, plabels = plabels)
      ## draw points
      if(any(object@adeg.par$ppoints$cex > 0))
        panel.points(x = c(xpoints, xpoints2), y = c(score1, score2), pch = object@adeg.par$ppoints$pch, cex = object@adeg.par$ppoints$cex, col = object@adeg.par$ppoints$col, alpha = object@adeg.par$ppoints$alpha, fill = object@adeg.par$ppoints$fill)
    }
  })


s1d.match <- function(score1, score2, labels = 1:NROW(score1), at = 0.5, facets = NULL, plot = TRUE, storeData = TRUE, add = FALSE, pos = -1, ...) {
  
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
      object <- multi.facets.S1(thecall, sortparameters$adepar, samelimits = sortparameters$g.args$samelimits)
    else 
      stop("Facets are not allowed with multiple scores")
  }
  
  ## multiple scores
  else if(NCOL(score1) > 1) { 
    object <- multi.score.S1(thecall)
  }
  
  ## simple ADEg graphic
  else {
    if(length(sortparameters$rest))
      warning(c("Unused parameters: ", paste(unique(names(sortparameters$rest)), " ", sep = "")), call. = FALSE)
    
    ## creation of the ADEg object
    if(storeData)
      tmp_data <- list(score = c(score1, score2), labels = labels, at = at, frame = sys.nframe() + pos, storeData = storeData)
    else
      tmp_data <- list(score = call("c", thecall$score1, thecall$score2), labels = thecall$labels, at = thecall$at, frame = sys.nframe() + pos, storeData = storeData)
    object <- new(Class = "S1.match", data = tmp_data, adeg.par = sortparameters$adepar, trellis.par = sortparameters$trellis, g.args = sortparameters$g.args, Call = match.call())
    
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
