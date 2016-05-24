###########################################################
##                          s1d.distri                   ##
###########################################################

setClass(
  Class = "S1.distri",
  contains = "ADEg.S1"
)


setMethod(
  f = "initialize",
  signature = "S1.distri",
  definition = function(.Object, data = list(score = NULL, dfdistri = NULL, labels = NULL, at = NULL, frame = 0, storeData = TRUE), ...) {
    .Object <- callNextMethod(.Object, data = data, ...) ## ADEg.S1 initialize
    .Object@data$dfdistri <- data$dfdistri
    .Object@data$labels <- data$labels
    return(.Object)
  })


setMethod(
  f = "prepare",
  signature = "S1.distri",
  definition = function(object) {
    name_obj <- deparse(substitute(object))
    
    ## pre-management of graphics parameters      
    oldparamadeg <- adegpar()
    on.exit(adegpar(oldparamadeg))
    adegtot <- adegpar(object@adeg.par)
    
    if(object@data$storeData) {
      dfdistri <- object@data$dfdistri
      score <- object@data$score
      labels <- object@data$labels      
  	} else {
      dfdistri <- eval(object@data$dfdistri, envir = sys.frame(object@data$frame))
      score <- eval(object@data$score, envir = sys.frame(object@data$frame))
      labels <- eval(object@data$labels, envir = sys.frame(object@data$frame))
    }
    
    score <- as.matrix(score)[, 1]  ## to manage 'score' when it is a data.frame with only one column
    
    ## change default for some parameters
    if(adegtot$p1d$horizontal & is.null(object@adeg.par$plabels$srt))
      adegtot$plabels$srt <- 0
    else if(!adegtot$p1d$horizontal & is.null(object@adeg.par$plabels$srt))
      adegtot$plabels$srt <- 90
    
    ## statistics calculus
    object@stats$means <- sapply(dfdistri, function(x) weighted.mean(score, x))
    names(object@stats$means) <- labels
    object@stats$sds <- sapply(dfdistri, function(x) sqrt(varwt(score, x)))
    
    ## object modification before calling inherited method
    object@adeg.par <- adegtot
    callNextMethod() ## prepare graph
    
    assign(name_obj, object, envir = parent.frame())
  })


setMethod(
  f = "panel",
  signature = "S1.distri",
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
    ngroups <- length(object@stats$means)
    means <- object@stats$means
    sds <- object@stats$sds * object@g.args$sdSize
    plabels <- object@adeg.par$plabels
    
    lead <- ifelse(pscore$reverse, -1, 1)

    if(pscore$horizontal) {
      ## horizontal plot
      ylab <- at
      if(object@g.args$yrank) {
        idx <- order(means, decreasing = TRUE)
        means <- means[idx]
        sds <- sds[idx]
        labels <- labels[idx]
      }
      
      do.call("panel.segments", c(list(x0 = means - sds, y0 = ylab, x1 = means + sds, y1 = ylab), object@adeg.par$plines))
      do.call("panel.points", c(list(x = means, y = ylab), object@adeg.par$ppoints))
      etis <- ifelse(abs(lims$xlim[2] - (means + sds)) > abs(lims$xlim[1] - (means - sds)), 1, -1)
    } else {
      ## vertical plot
      xlab <- at
      if(object@g.args$yrank) {
        idx <- order(means, decreasing = TRUE)
        means <- means[idx]
        sds <- sds[idx]
        labels <- labels[idx]
      }

      do.call("panel.segments", c(list(x0 = xlab, y0 = means - sds, x1 = xlab, y1 = means + sds), object@adeg.par$plines))
      do.call("panel.points", c(list(x = xlab, y = means), object@adeg.par$ppoints))
      etis <- ifelse(abs(lims$ylim[2] - (means + sds)) > abs(lims$ylim[1] - (means - sds)), 1, -1)
    }
    
    ## manage string rotation
    srt <- 0
    if(is.numeric(plabels$srt[1]))
      srt <- plabels$srt[1]
    else {
      if(plabels$srt[1] == "horizontal")
        srt <- 0
      else if(plabels$srt[1] == "vertical")
        srt <- 90
    }
    
    ## draw labels
    if(abs(sin(srt)) > sin(45)) {
      ## almost vertical labels
      if(pscore$horizontal)
        width <- stringWidth("h")
      else
        width <- stringWidth(labels) + stringWidth("h")
      
      width <- rep(plabels$cex, length.out = length(labels)) * convertUnit(width, "native", typeFrom = "dimension", axisFrom = "x", axisTo = "y", valueOnly = TRUE) / 2 
    } else {
      ## almost horizont labels
      if(pscore$horizontal)
        width <- stringWidth(labels) + stringWidth("h")
      else
        width <- stringWidth("h")
      
      width <- rep(plabels$cex, length.out = length(labels)) * convertUnit(width, "native", typeFrom = "dimension", axisFrom = "x", valueOnly = TRUE) / 2 
    }

    if(pscore$horizontal)
      adeg.panel.label(x = means + etis * (sds + width), y = ylab, labels = labels, plabels = plabels)
    else
      adeg.panel.label(x = xlab, y = means + etis * (sds + width), labels = labels, plabels = plabels)
  })


s1d.distri <- function(score, dfdistri, labels = colnames(dfdistri), at = 1:NCOL(dfdistri), yrank = TRUE, sdSize = 1, facets = NULL, plot = TRUE, storeData = TRUE, add = FALSE, pos = -1, ...) {

  ## evaluation of some parameters (required for multiplot)
  thecall <- .expand.call(match.call())
  dfdistri <- eval(thecall$dfdistri, envir = sys.frame(sys.nframe() + pos))
  score <- eval(thecall$score, envir = sys.frame(sys.nframe() + pos))
  if(NROW(dfdistri) != NROW(score))
    stop("dfdis and score must have the same number of rows")
  
  ## parameters sorted
  sortparameters <- sortparamADEg(...)
  
  ## facets
  if(!is.null(facets)){
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
    g.args <- c(sortparameters$g.args, list(yrank = yrank, sdSize = sdSize))
    if(storeData)
    	tmp_data <- list(score = score, dfdistri = dfdistri, at = at, labels = labels, frame = sys.nframe() + pos, storeData = storeData)
    else
      tmp_data <- list(score = thecall$score, dfdistri = thecall$dfdistri, at = thecall$at, labels = labels, frame = sys.nframe() + pos, storeData = storeData)
    object <- new(Class = "S1.distri", data = tmp_data, adeg.par = sortparameters$adepar, trellis.par = sortparameters$trellis, g.args = g.args, Call = match.call())
    
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
