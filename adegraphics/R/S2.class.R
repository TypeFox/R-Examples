##########################################################################
##                            s.class                                   ##
##########################################################################

setClass(
  Class = "S2.class",
  contains = "ADEg.S2",
)

  
setMethod(
  f = "initialize",
  signature = "S2.class",
  definition = function(.Object, data = list(dfxy = NULL, xax = 1, yax = 2, fac = NULL, wt = NULL, labels = NULL, frame = 0, storeData = TRUE),  ...) {
    .Object <- callNextMethod(.Object, data = data, ...)
    .Object@data$fac <- data$fac
    .Object@data$wt <- data$wt
    .Object@data$labels <- data$labels
    return(.Object)
  })


setMethod(
  ## prepare computations for ellipses, stars and labels
  f = "prepare",
  signature = "S2.class",
  definition = function(object) {
    name_obj <- deparse(substitute(object))
    
    ## pre-management of graphics parameters      
    oldparamadeg <- adegpar()
    on.exit(adegpar(oldparamadeg))
    adegtot <- adegpar(object@adeg.par)
    
    if(object@data$storeData) {
      fac <- as.factor(object@data$fac)
      dfxy <- object@data$dfxy
      wt <- object@data$wt
    } else {
      fac <- as.factor(eval(object@data$fac, envir = sys.frame(object@data$frame)))
      dfxy <- eval(object@data$dfxy, envir = sys.frame(object@data$frame))
      wt <- eval(object@data$wt, envir = sys.frame(object@data$frame))
    }
    
    if(is.null(object@adeg.par$porigin$include) & (any(names(object@g.args) %in% c("Sp", "nbobject"))))
      adegtot$porigin$include <- FALSE

    if(any(adegtot$plabels$cex > 0) & is.null(object@adeg.par$plegend$drawKey)) ## if labels, no legend
        adegtot$plegend$drawKey <- FALSE
    ## setting colors 
    if(!is.null(object@g.args$col)){
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
        
        if(col.idx){
            if(is.null(object@adeg.par$ppoints$col))
                adegtot$ppoints$col <- colT
            if(is.null(object@adeg.par$ppoints$fill))
                adegtot$ppoints$fill <- colT
            if(is.null(object@adeg.par$pellipses$border))
                adegtot$pellipses$border <- colT
            if(is.null(object@adeg.par$pellipses$col))
                adegtot$pellipses$col <- colT
            if(is.null(object@adeg.par$plabels$col))
                adegtot$plabels$col <- colT
            if(is.null(object@adeg.par$plabels$boxes$border))
                adegtot$plabels$boxes$border <- colT
            if(is.null(object@adeg.par$ppolygons$border))
                adegtot$ppolygons$border <- colT
            if(is.null(object@adeg.par$ppolygons$col))
                adegtot$ppolygons$col <- colT
            if(is.null(object@adeg.par$plines$col))
                adegtot$plines$col <- colT
        }
    }
    
    ## preliminary computations
    object@stats$means <- matrix(meanfacwt(dfxy[, c(object@data$xax, object@data$yax)], fac, wt), nrow = nlevels(fac))
    ## for ellipse, covariance and variance needed
    if(object@g.args$ellipseSize)
      object@stats$covvar <- covfacwt(dfxy[, c(object@data$xax, object@data$yax)], fac, wt)
    
    ## object modification before calling inherited method
    object@adeg.par <- adegtot
    callNextMethod() ## prepare graph
    
    ## compute ellipses
    if(object@g.args$ellipseSize > 0) { 
      object@s.misc$ellipses <- lapply(1:nlevels(fac), FUN = function(i) {
        .util.ellipse(object@stats$means[i, 1], object@stats$means[i, 2], vx = object@stats$covvar[[i]][1, 1], vy = object@stats$covvar[[i]][2, 2],
                      cxy = object@stats$covvar[[i]][1, 2], coeff = object@g.args$ellipseSize)
      })                                                    
    }
    
    ## compute convex hulls
    if(!is.null(object@g.args$chullSize))
      if(any(object@g.args$chullSize > 0))
        object@s.misc$chullcoord <- .util.chull(dfxy[, object@data$xax], dfxy[, object@data$yax], object@stats$means[, 1], object@stats$means[, 2], fac = fac, chullSize = object@g.args$chullSize)
    
    ## never optimized labels for s.class
    object@adeg.par$plabels$optim <- FALSE
    assign(name_obj, object, envir = parent.frame())
  })


## a changer: dessin level par level, 
setMethod(
  f = "panel",
  signature = "S2.class",
  definition = function(object, x, y) {
    if(object@data$storeData) {
      fac <- object@data$fac
      labels <- object@data$labels    
    } else {
      fac <- eval(object@data$fac, envir = sys.frame(object@data$frame))
      labels <- eval(object@data$labels, envir = sys.frame(object@data$frame))
    }
    nlev <- nlevels(fac)
    
    ## convex hulls
    if(any(object@g.args$chullSize > 0)) {
      chullpo <- object@s.misc$chullcoord
      ppolygons <- lapply(object@adeg.par$ppolygons, FUN = function(x) {rep(x, length.out = length(chullpo))})
      
      for(level in 1:nlev) 
        if(!any(is.null(chullpo[[level]]))) {
          for(j in 1:length(chullpo[[level]]))
            panel.polygon(
              x = chullpo[[level]][[j]][, 1], y = chullpo[[level]][[j]][, 2],
              border = ppolygons$border[level], col = ppolygons$col[level],
              lty = ppolygons$lty[level], lwd = ppolygons$lwd[level], alpha = ppolygons$alpha[level])
        }
    }
    ## ellipses
    if(object@g.args$ellipseSize > 0) {
      ellip <- object@s.misc$ellipses
      pellip <- object@adeg.par$pellipses
      pellip <- lapply(pellip, FUN = function(x) {if(is.list(x)) return(x) else rep(x, le = length(ellip))})
      pellip$axes <- lapply(pellip$axes, FUN = function(x) {rep(x, length.out = length(ellip))})
      
      for(level in 1:nlev) { ## for each group
        ell <- ellip[[level]]
        if(!(any(is.null(ell))))
          if(!any(is.na(ell))) {
            panel.polygon(ell$x, ell$y, col = pellip$col[level], lwd = pellip$lwd[level], lty = pellip$lty[level], alpha = pellip$alpha[level], border = pellip$border[level])
            if(pellip$axes$draw[level]) { ## axes drawing
              panel.segments(ell$seg1[1], ell$seg1[2], ell$seg1[3], ell$seg1[4], lwd = pellip$axes$lwd[level], lty = pellip$axes$lty[level], col = pellip$axes$col[level])
              panel.segments(ell$seg2[1], ell$seg2[2], ell$seg2[3], ell$seg2[4], lwd = pellip$axes$lwd[level], lty = pellip$axes$lty[level], col = pellip$axes$col[level])
            }
          }
      }
    }
    
    ## stars
    if(object@g.args$starSize > 0) {
      plines <- lapply(object@adeg.par$plines, FUN = function(x) {rep(x, length.out = nlev)})
      for(level in 1:nlev) {
        if(all(is.finite(object@stats$means[level, ]))) {
          xbase <- object@stats$means[level, 1]
          ybase <- object@stats$means[level, 2]
          xlev <- x[fac == levels(fac)[level]]
          ylev <- y[fac == levels(fac)[level]]
          panel.segments(
            x0 = xbase, 
            y0 = ybase,
            x1 = xbase + object@g.args$starSize * (xlev - xbase),
            y1 = ybase + object@g.args$starSize * (ylev - ybase),
            lty = plines$lty[level], lwd = plines$lwd[level], col = plines$col[level])
        }
      }
    }
    
    ## plot points
    if(any(object@adeg.par$ppoints$cex > 0)) {
      ppoints <- object@adeg.par$ppoints
      if(nlev > 1) {
        ppoints <- lapply(object@adeg.par$ppoints, FUN = function(x, fac, nlev) {
          if(length(x) > nlev)
            return(x)
          else {
            xlev <- rep(x, length.out = nlev)
            xpar <- xlev[fac]
            return(xpar)
          }
        }, fac = fac, nlev = nlev)
      }
      
      if(any(is.na(ppoints$pch))) {
        indx <- 1:length(x)
        indx <- indx[- which(is.na(ppoints$pch))]
        panel.points(x = x[indx], y = y[indx], type = "p", pch = ppoints$pch[indx], cex = ppoints$cex[indx],
                     col = ppoints$col[indx], alpha = ppoints$alpha[indx], fill = ppoints$fill[indx])
      } else
        panel.points(x = x, y = y, type = "p", pch = ppoints$pch, cex = ppoints$cex, col = ppoints$col,
                     alpha = ppoints$alpha, fill = ppoints$fill)
    }
    
    ## plot of labels
    if(any(object@adeg.par$plabels$cex > 0)) {
      labX <- object@stats$means[, 1]
      labY <- object@stats$means[, 2]
      adeg.panel.label(x = labX, y = labY, labels = labels, object@adeg.par$plabels)
    }
  })

  
s.class <- function(dfxy, fac, xax = 1, yax = 2, wt = rep(1, NROW(fac)), labels = levels(fac), ellipseSize = 1.5, starSize = 1, 
  									chullSize = NULL, col = NULL, facets = NULL, plot = TRUE, storeData = TRUE, add = FALSE, pos = -1, ...) {
  
  ## evaluation of some parameters (required for multiplot)
  thecall <- .expand.call(match.call())
  labels <- eval(thecall$labels, envir = sys.frame(sys.nframe() + pos))
  fac <- eval(thecall$fac, envir = sys.frame(sys.nframe() + pos))
  df <- try(as.data.frame(eval(thecall$dfxy, envir = sys.frame(sys.nframe() + pos))), silent = TRUE)
  if((class(df) == "try-error") | is.null(thecall$dfxy)) ## non convenient dfxy argument
    stop("non convenient selection for dfxy (can not be converted to dataframe)")
  if(missing(fac))
    stop("no factor specified")
  
  if(NCOL(fac) == 1) {
    fac <- as.factor(fac)
    if(length(labels) != nlevels(fac))
      stop("wrong number of labels")
  }
  
  ## parameters sorted
  sortparameters <- sortparamADEg(...)
  
  ## facets
  if(!is.null(facets)) { 
    if((length(xax) == 1 & length(yax) == 1) & NCOL(fac) == 1)
      object <- multi.facets.S2(thecall, sortparameters$adepar, samelimits = sortparameters$g.args$samelimits)
    else 
      stop("Facets are not allowed with multiple xax/yax or multiple fac")
  }
  
  ## multiple axes
  else if((length(xax) > 1 | length(yax) > 1)) {
    if(NCOL(fac) == 1)
      object <- multi.ax.S2(thecall)
    else 
      stop("Multiple xax/yax are not allowed with multiple fac")
  }
  
  ## multiple fac
  else if(NCOL(fac) > 1) {
    object <- multi.variables.S2(thecall, "fac")
  }
  
  ## simple ADEg graphic
  else {
    if(length(sortparameters$rest))
  	  warning(c("Unused parameters: ", paste(unique(names(sortparameters$rest)), " ", sep = "")), call. = FALSE)
    
    ## creation of the ADEg object
    g.args <- c(sortparameters$g.args, list(ellipseSize = ellipseSize, starSize = starSize, chullSize = chullSize, col = col))
    if(storeData)
    	tmp_data <- list(dfxy = dfxy, fac = fac, xax = xax, yax = yax, wt = wt, labels = labels, frame = sys.nframe() + pos, storeData = storeData)
    else
      tmp_data <- list(dfxy = thecall$dfxy, fac = thecall$fac, xax = xax, yax = yax, wt = thecall$wt, labels = thecall$labels, frame = sys.nframe() + pos, storeData = storeData)
    object <- new(Class = "S2.class", data = tmp_data, adeg.par = sortparameters$adepar, trellis.par = sortparameters$trellis, g.args = g.args, Call = as.call(thecall))

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
