###########################################################
##                          s1d.class                    ##
###########################################################

setClass(
  Class = "S1.class",
  contains = "ADEg.S1",
)


setMethod(
  f = "initialize",
  signature = "S1.class",
  definition = function(.Object, data = list(score = NULL, fac = NULL, wt = NULL, labels = NULL, at = NULL, frame = 0, storeData = TRUE), ...) {
    .Object <- callNextMethod(.Object, data = data, ...) ## ADEg.S1 initialize
    .Object@data$fac <- data$fac
    .Object@data$wt <- data$wt
    .Object@data$labels <- data$labels
    return(.Object)
  })


setMethod(
  f = "prepare",
  signature = "S1.class",
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
    
    ## change default for some parameters
    if(adegtot$p1d$horizontal & is.null(object@adeg.par$plabels$srt))
      adegtot$plabels$srt <- 90
    else if(!adegtot$p1d$horizontal & is.null(object@adeg.par$plabels$srt))
      adegtot$plabels$srt <- 0
    
    if(any(adegtot$plabels$cex > 0) & is.null(object@adeg.par$plegend$drawKey)) ## if labels, no legend
      adegtot$plegend$drawKey <- FALSE    
    ## setting colors 
    if(!is.null(object@g.args$col)) {
      col.idx <- FALSE
      if(is.logical(object@g.args$col)) {
        if(object@g.args$col){
          colT <- adegtot$ppalette$quali(nlevels(fac))
          col.idx <- TRUE
        }
      } else {
        colT <- rep(object@g.args$col, length.out = nlevels(fac))
        col.idx <- TRUE
      }
      
      if(col.idx) {
        if(is.null(object@adeg.par$ppoints$col))
          adegtot$ppoints$col <- colT
        if(is.null(object@adeg.par$ppoints$fill))
          adegtot$ppoints$fill <- colT
        if(is.null(object@adeg.par$plabels$col))
          adegtot$plabels$col <- colT
        if(is.null(object@adeg.par$plabels$boxes$border))
          adegtot$plabels$boxes$border <- colT
        if(is.null(object@adeg.par$plines$col))
          adegtot$plines$col <- colT
      }
    }
    
    if(adegtot$p1d$horizontal & is.null(object@g.args$ylim))
      object@g.args$ylim <- c(0, 1)
    
    if(!adegtot$p1d$horizontal & is.null(object@g.args$xlim))
      object@g.args$xlim <- c(0, 1)
    
    ## object modification before calling inherited method
    object@adeg.par <- adegtot
    callNextMethod() ## prepare graph
    
    assign(name_obj, object, envir = parent.frame())
  })


## TODO: label orientation (works only for horizontal / vertical labels)
setMethod(
  f= "panel",
  signature = "S1.class",
  definition = function(object, x, y) {
    
    if(object@data$storeData) {
      fac <- object@data$fac
      score <- object@data$score
      wt <- object@data$wt
      at <- object@data$at
      labels <- object@data$labels
    } else {
      fac <- eval(object@data$fac, envir = sys.frame(object@data$frame))
      score <- eval(object@data$score, envir = sys.frame(object@data$frame))
      wt <- eval(object@data$wt, envir = sys.frame(object@data$frame))
      at <- eval(object@data$at, envir = sys.frame(object@data$frame))
      labels <- eval(object@data$labels, envir = sys.frame(object@data$frame))
    }
    
    fac <- as.factor(fac)
    nlev <- nlevels(fac)
    object@stats$means <- meanfacwt(score, fac, wt = wt)
    lims <- current.panel.limits(unit = "native")
    pscore <- object@adeg.par$p1d
    ## repeat graphical parameters (one for each level)
    ppoints <- lapply(object@adeg.par$ppoints, FUN = function(x) x <- rep(x, length.out = nlev))
    plines <- lapply(object@adeg.par$plines, FUN = function(x) x <- rep(x, length.out = nlev))
    plabels <- lapply(object@adeg.par$plabels, FUN = function(x) x <- rep(x, length.out = nlev))
    plboxes <- lapply(object@adeg.par$plabels$boxes, FUN = function(x) x <- rep(x, length.out = nlev))
    plabels$boxes <- plboxes
    
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
      
      ## get positions for labels
      if(object@g.args$poslabel == "regular") {
        spacelab <- diff(lims$xlim) / (nlev + 1)
        xlab <- seq(from = lims$xlim[1] + spacelab, by = spacelab, length.out = nlev)[rank(object@stats$means, ties.method = "first")]
      } else
        xlab <- object@stats$means
      
      ## repeat means for each individual   
      xlablines <- xlab[fac]
      
      ## repeat ylab for each individual
      ylab <- rep(at, length.out = nlev)
      ylablines <- ylab[fac]
      
      ## draw lines and labels
      ypoints <- object@s.misc$rug
      panel.segments(x0 = xpoints, y0 = ypoints, x1 = xlablines, y1 = ylablines, lwd = plines$lwd[fac], col = plines$col[fac], lty = plines$lty[fac])
      if(any(ppoints$cex > 0))
        panel.points(x = xpoints, y = ypoints, pch = ppoints$pch[fac], cex = ppoints$cex[fac], col = ppoints$col[fac], alpha = ppoints$alpha[fac], fill = ppoints$fill[fac])
      if(any(plabels$cex > 0))
        adeg.panel.label(x = xlab, y = ylab + lead * h / 2, labels = labels, plabels = plabels)
      
    } else {
      ## vertical plot
      ypoints <- y
      
      ## get positions for labels
      if(object@g.args$poslabel == "regular") {
        spacelab <- diff(lims$ylim) / (nlev + 1)
        ylab <- seq(from = lims$ylim[1] + spacelab, by = spacelab, length.out = nlev)[rank(object@stats$means, ties.method = "first")]
      } else
        ylab <- object@stats$means
      
      ## repeat means for each individual   
      ylablines <- ylab[fac]
      
      ## repeat ylab for each individual
      xlab <- rep(at, length.out = nlev)
      xlablines <- xlab[fac]
      
      ## draw lines and labels
      xpoints <- object@s.misc$rug
      panel.segments(x0 = xpoints, y0 = ypoints, x1 = xlablines, y1 = ylablines, lwd = plines$lwd[fac], col = plines$col[fac], lty = plines$lty[fac])
      if(any(ppoints$cex > 0))
        panel.points(x = xpoints, y = ypoints, pch = ppoints$pch[fac], cex = ppoints$cex[fac], col = ppoints$col[fac], alpha = ppoints$alpha[fac], fill = ppoints$fill[fac])
      if(any(plabels$cex > 0))
        adeg.panel.label(x = xlab + lead * w / 2 , y = ylab, labels = labels, plabels = plabels)
    } 
  })


s1d.class <- function(score, fac, wt = rep(1, NROW(fac)), labels = levels(fac), at = 0.5, poslabel = c("regular", "value"), col = NULL, 
                      facets = NULL, plot = TRUE, storeData = TRUE, add = FALSE, pos = -1, ...) {
  
  ## evaluation of some parameters
  thecall <- .expand.call(match.call())
  labels <- eval(thecall$labels, envir = sys.frame(sys.nframe() + pos))
  fac <- eval(thecall$fac, envir = sys.frame(sys.nframe() + pos))
  score <- eval(thecall$score, envir = sys.frame(sys.nframe() + pos))
  if(NCOL(fac) == 1) {
    fac <- as.factor(fac)
    if(length(labels) != nlevels(fac))
      stop("wrong number of labels")
  }
  if(NROW(score) != NROW(fac))
    stop("score and factor must have the same number of rows")
  
  ## parameters sorted
  sortparameters <- sortparamADEg(...)
  
  ## facets
  if(!is.null(facets)) {
    if(NCOL(score) == 1 & NCOL(fac) == 1)
      object <- multi.facets.S1(thecall, sortparameters$adepar, samelimits = sortparameters$g.args$samelimits)
    else 
      stop("Facets are not allowed with multiple scores or fac")
  }
  
  ## multiple scores
  else if(NCOL(score) > 1) {
    if(NCOL(fac) == 1)
      object <- multi.score.S1(thecall)
    else 
      stop("Multiple scores are not allowed with multiple fac")
  }
  
  ## multiple fac
  else if(NCOL(fac) > 1) {
    object <- multi.variables.S1(thecall, "fac")
  }
  
  ## simple ADEg graphic
  else {
    if(length(sortparameters$rest))
      warning(c("Unused parameters: ", paste(unique(names(sortparameters$rest)), " ", sep = "")), call. = FALSE)
    
    ## creation of the ADEg object
    g.args <- c(sortparameters$g.args, list(poslabel = match.arg(poslabel), col = col))
    if(storeData)
      tmp_data <- list(score = score, wt = wt, fac = fac, labels = labels, at = at, frame = sys.nframe() + pos, storeData = storeData)
    else
      tmp_data <- list(score = thecall$score, wt = thecall$wt, fac = thecall$fac, labels = thecall$labels, at = thecall$at, frame = sys.nframe() + pos, storeData = storeData)
    object <- new(Class = "S1.class", data = tmp_data, adeg.par = sortparameters$adepar, trellis.par = sortparameters$trellis, g.args = g.args, Call = match.call())
    
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

