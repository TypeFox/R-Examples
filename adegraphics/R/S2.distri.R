##########################################################################
##                           s.distri                                   ##
##########################################################################

setClass(
  Class = "S2.distri",
  contains = "ADEg.S2",
)


setMethod(
  f = "initialize",
  signature = "S2.distri",
  definition = function(.Object, data = list(dfxy = NULL, dfdistri = NULL, xax = 1, yax = 2, frame = 0, storeData = TRUE), ...) {
    .Object <- callNextMethod(.Object, data = data, ...)
    .Object@data$dfdistri <- data$dfdistri
    return(.Object)
  })


setMethod(
  ## prepare computations for ellipses, stars and labels
  f = "prepare",
  signature = "S2.distri",
  definition = function(object) {
    name_obj <- deparse(substitute(object))
    
    ## pre-management of graphics parameters      
    oldparamadeg <- adegpar()
    on.exit(adegpar(oldparamadeg))
    adegtot <- adegpar(object@adeg.par)
    
    if(!object@data$storeData) {
      dfxy <- eval(object@data$dfxy, envir = sys.frame(object@data$frame))
      dfdistri <- eval(object@data$dfdistri, envir = sys.frame(object@data$frame))
  	} else {
      dfxy <- object@data$dfxy
      dfdistri <- object@data$dfdistri
    }
    
    if(is.null(colnames(dfdistri))) 
      adegtot$plabels$cex <- 0 ## no labels if no colnames in original data
    if(is.null(object@adeg.par$porigin$include) & (any(names(object@g.args) %in% c("Sp", "nbobject"))))
      adegtot$porigin$include <- FALSE
    
    ## setting colors
    if(!is.null(object@g.args$col))
      if(is.logical(object@g.args$col)) {
        if(object@g.args$col)
          adegtot$pellipses$border <- adegtot$pellipses$col <- adegtot$plabels$col <- adegtot$plabels$boxes$border <- adegtot$plines$col <- adegtot$ppalette$quali(NCOL(dfdistri))
      } else
        adegtot$pellipses$border <- adegtot$pellipses$col <- adegtot$plabels$col <- adegtot$plabels$boxes$border <- adegtot$plines$col <- rep(object@g.args$col, length.out = NCOL(dfdistri))
    
    ## statistics calculus
    object@stats$means <- t(apply(as.data.frame(dfdistri), 2, FUN = function(x) {apply(dfxy[, c(object@data$xax, object@data$yax)], 2, weighted.mean , w = x)}))
    if(object@g.args$ellipseSize)
      object@stats$covvar <- lapply(as.data.frame(dfdistri), FUN = function(x) {covwt(dfxy[, c(object@data$xax, object@data$yax)], wt = x)})
    else 
      object@stats$covvar <- NULL
    
    ## object modification before calling inherited method
    object@adeg.par <- adegtot
    callNextMethod() ## prepare graph
    
    ## compute ellipses
    if(object@g.args$ellipseSize > 0) { 
      object@s.misc$ellipses <- lapply(1:nrow(object@stats$means), FUN = function(i) {
        .util.ellipse(object@stats$means[i, 1], object@stats$means[i, 2], vx = object@stats$covvar[[i]][1, 1], vy = object@stats$covvar[[i]][2, 2],
                      cxy = object@stats$covvar[[i]][1, 2], coeff = object@g.args$ellipseSize)
      })                                                    
    }
    ## never optimized labels for s.distri
    object@adeg.par$plabels$optim <- FALSE
    assign(name_obj, object, envir = parent.frame())
  })


setMethod(
  f = "panel",
  signature = "S2.distri",
  definition = function(object, x, y) {
    
    if(object@data$storeData)
      dfdistri <- object@data$dfdistri
    else
      dfdistri <- eval(object@data$dfdistri, envir = sys.frame(object@data$frame))
    
    ## ellipses
    if(object@g.args$ellipseSize > 0) {
      ellip <- object@s.misc$ellipses
      pellip <- object@adeg.par$pellipses
      pellip <- lapply(pellip, FUN = function(x) {if(is.list(x)) return(x) else rep(x, length.out = length(ellip))})
      pellip$axes <- lapply(pellip$axes, FUN = function(x) {rep(x, length.out = length(ellip))})
      
      for(group in 1:NCOL(dfdistri)) { ## for each group
        ell <- ellip[[group]]
        if(!(any(is.null(ell))))
          if(!any(is.na(ell))) {
            panel.polygon(ell$x, ell$y, col = pellip$col[group], lwd = pellip$lwd[group],
                          lty = pellip$lty[group], alpha = pellip$alpha[group], border = pellip$border[group])
            if(pellip$axes$draw[group]) {
              ## axes drawing
              panel.segments(ell$seg1[1], ell$seg1[2], ell$seg1[3], ell$seg1[4], lwd = pellip$axes$lwd[group],
                             lty = pellip$axes$lty[group], col = pellip$axes$col[group])
              panel.segments(ell$seg2[1], ell$seg2[2], ell$seg2[3], ell$seg2[4], lwd = pellip$axes$lwd[group],
                             lty = pellip$axes$lty[group], col = pellip$axes$col[group])
            }
          }
      }
    }
    
    ## stars
    if(object@g.args$starSize > 0) {
      plines <- lapply(object@adeg.par$plines, FUN = function(x) {rep(x, length.out = NCOL(dfdistri))})
      for(group in 1:NCOL(dfdistri)) {
        if(all(is.finite(object@stats$means[group, ]))) {
          xbase <- object@stats$means[group, 1]
          ybase <- object@stats$means[group, 2]
          xlev <- x[which(as.data.frame(dfdistri)[, group] > 0)]
          ylev <- y[which(as.data.frame(dfdistri)[, group] > 0)]
          panel.segments(
            x0 = xbase,
            y0 = ybase,
            x1 = xbase + object@g.args$starSize * (xlev - xbase),
            y1 = ybase + object@g.args$starSize * (ylev - ybase),
            lty = plines$lty[group], lwd = plines$lwd[group], col = plines$col[group])
        } 
      }
    }
    
    ## plot points
    if(any(object@adeg.par$ppoints$cex > 0)) {
      ppoints <- lapply(object@adeg.par$ppoints, function(x) rep(x, length.out = NROW(dfdistri)))
      
      if(any(is.na(ppoints$pch))) {
        indx <- 1:length(x)
        indx <- indx[- which(is.na(ppoints$pch))]
        panel.points(x = x[indx], y = y[indx], type = "p", pch = ppoints$pch[indx], cex = ppoints$cex[indx],
                     col = ppoints$col[indx], alpha = ppoints$alpha[indx], fill = ppoints$fill[indx])}
      else
        panel.points(x = x, y = y, type = "p", pch = ppoints$pch, cex = ppoints$cex, col = ppoints$col,
                     alpha = ppoints$alpha, fill = ppoints$fill)
    }
    ## plot of labels
    if(any(object@adeg.par$plabels$cex > 0)) {
      labX <- object@stats$means[, 1]
      labY <- object@stats$means[, 2]
      adeg.panel.label(x = labX, y = labY, labels = colnames(dfdistri), object@adeg.par$plabels)
    }
  })
  

s.distri <- function(dfxy, dfdistri, xax = 1, yax = 2, starSize = 1, ellipseSize = 1.5, col = NULL, facets = NULL, plot = TRUE, 
  storeData = TRUE, add = FALSE, pos = -1, ...) {
  
  ## evaluation of some parameters (required for multiplot)
  thecall <- .expand.call(match.call())
  df <- try(as.data.frame(eval(thecall$dfxy, envir = sys.frame(sys.nframe() + pos))), silent = TRUE)
  if((class(df) == "try-error") | is.null(dfxy))
    stop("dfxy, can not be converted as dataframe or is NULL")
  
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
    g.args <- c(sortparameters$g.args, list(ellipseSize = ellipseSize, starSize = starSize, col = col))
    if(storeData)
      tmp_data <- list(dfxy = dfxy, dfdistri = dfdistri, xax = xax, yax = yax, frame = sys.nframe() + pos, storeData = storeData)
    else
    	tmp_data <- list(dfxy = thecall$dfxy, dfdistri = thecall$dfdistri, xax = xax, yax = yax, frame = sys.nframe() + pos, storeData = storeData)
    object <- new(Class = "S2.distri", data = tmp_data, adeg.par = sortparameters$adepar, trellis.par = sortparameters$trellis, g.args = g.args, Call = as.call(thecall))

    ## preparation of the graph
    prepare(object)
    setlatticecall(object)
    if(add)
      object <- add.ADEg(object)
  }
  
  if(! add & plot)
    print(object)
  invisible(object)
}

