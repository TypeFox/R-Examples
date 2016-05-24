######################################################
##                   Tr.class                      ###
##    Triangular representation with a factor      ###
######################################################

setClass(
  Class = "Tr.class",
  contains = "ADEg.Tr"
)


setMethod(
  f = "initialize",
  signature = "Tr.class",
  definition = function(.Object, data = list(dfxyz = NULL, fac = NULL, wt = NULL, labels = NULL, frame = 0, storeData = TRUE), ...) {
    .Object <- callNextMethod(.Object, data = data, ...)
    .Object@data$fac <- data$fac
    .Object@data$wt <- data$wt
    .Object@data$labels <- data$labels
    return(.Object)
  })


setMethod(
  f = "prepare",
  signature = "Tr.class",
  definition = function(object) {
    name_obj <- deparse(substitute(object))
    
    if(object@data$storeData) {
      df <- object@data$dfxyz
      fac <- as.factor(object@data$fac)
      wt <- object@data$wt
    } else {
      fac <- as.factor(eval(object@data$fac, envir = sys.frame(object@data$frame)))
      df <- eval(object@data$dfxyz, envir = sys.frame(object@data$frame))
      wt <- eval(object@data$wt, envir = sys.frame(object@data$frame))
    }
    nlev <- nlevels(fac)
    
    ## pre-management of graphics parameters      
    oldparamadeg <- adegpar()
    on.exit(adegpar(oldparamadeg))
    adegtot <- adegpar(object@adeg.par)
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
    
    ## object modification before calling inherited method
    object@adeg.par <- adegtot
    callNextMethod() ## prepare graph
    
    ## calculate 2D coordinates
    df <- sweep(df, 1, rowSums(df), "/")
    object@stats$coords2d <- .coordtotriangleM(df, mini3 = object@g.args$min3d, maxi3 = object@g.args$max3d)[, 2:3]
    
    ## compute means for the 3 variables (for getstats)
    object@stats$means <- matrix(meanfacwt(df, fac, wt), nrow = nlev)
    ## mean2d: columns: axes, row: levels
    object@stats$mean2d <- matrix(meanfacwt(object@stats$coords2d, fac, wt), nrow = nlev)
    mean.x <- object@stats$mean2d[, 1] ## all means rows as levels, columns as variables
    mean.y <- object@stats$mean2d[, 2]
    
    ## ellipses
    if(object@g.args$ellipseSize > 0) {
      object@stats$covvar <- covfacwt(df, fac, wt)
      object@stats$covvar2d <- covfacwt(object@stats$coords2d, fac, wt)
      covvartotal <- object@stats$covvar2d
      
      object@s.misc$ellipses <- lapply(1:nlev,
                                       FUN = function(i) {
                                         .util.ellipse(mean.x[i], mean.y[i], vx = covvartotal[[i]][1, 1], vy = covvartotal[[i]][2, 2], cxy = covvartotal[[i]][1, 2],
                                                       coeff = object@g.args$ellipseSize)
                                       })
    }

    ## convex hull
    if(!is.null(object@g.args$chullSize))
      if(any(object@g.args$chullSize > 0))
        object@s.misc$chullcoord  <- .util.chull(object@stats$coords2d[, 1], object@stats$coords2d[, 2], mean.x, mean.y, fac = fac, chullSize =  object@g.args$chullSize)
    
    object@adeg.par$plabels$optim <- FALSE
    assign(name_obj, object, envir = parent.frame())
  })


setMethod(
  f = "panel",
  signature = "Tr.class",
  definition = function(object, x, y) {

    if(object@data$storeData) {
      df <- object@data$dfxyz
      fac <- object@data$fac
      labels <- object@data$labels
    } else {
      fac <- eval(object@data$fac, envir = sys.frame(object@data$frame))
      df <- eval(object@data$dfxyz, envir = sys.frame(object@data$frame))
      labels <- eval(object@data$labels, envir = sys.frame(object@data$frame))
    }

    fac <- as.factor(fac)
    nlev <- nlevels(fac)
    
    ## draw convex hulls
    if(any(object@g.args$chullSize > 0)) {
      chullpo <- object@s.misc$chullcoord
      ppolygons <- lapply(object@adeg.par$ppolygons, FUN = function(x) rep(x, length.out = length(chullpo)))
      for(level in 1:nlev) {
        chull <- chullpo[[level]]
        for(j in 1:length(chull))
          panel.polygon(x = chull[[j]][, 1], y = chull[[j]][, 2], border = ppolygons$border[level], col = ppolygons$col[level], lty = ppolygons$lty[level], lwd = ppolygons$lwd[level], alpha = ppolygons$alpha[level])
      }}

    ## draw ellipses
    if(object@g.args$ellipseSize > 0) {
      ellip <- object@s.misc$ellipses
      pellip <- object@adeg.par$pellipses
      ## setting parameters, number of levels
      pellip <- lapply(pellip, FUN = function(x) {if(is.list(x)) return(x) else rep(x, length.out = length(ellip))})
      pellip$axes <- lapply(pellip$axes, FUN = function(x) {rep(x, length.out = length(ellip))})
      for(level in 1:nlev) {
        ell <- ellip[[level]]
        if(!(any(is.null(ell))))
          if(!any(is.na(ell))) {
            panel.polygon(ell$x, ell$y, col = pellip$col[level], lwd = pellip$lwd[level], lty = pellip$lty[level], alpha = pellip$alpha[level], border = pellip$border[level])
            if(pellip$axes$draw[level]) {
            	## draw axes
              panel.segments(ell$seg1[1], ell$seg1[2], ell$seg1[3], ell$seg1[4], lwd = pellip$axes$lwd[level], lty = pellip$axes$lty[level], col = pellip$axes$col[level])
              panel.segments(ell$seg2[1], ell$seg2[2], ell$seg2[3], ell$seg2[4], lwd = pellip$axes$lwd[level], lty = pellip$axes$lty[level], col = pellip$axes$col[level])
            }
          }
      }
    }
    
    ## draw stars
    if(object@g.args$starSize > 0) {
      plines <- lapply(object@adeg.par$plines, FUN = function(x) {rep(x, length.out = nlev)})
      xlx <- split(object@stats$coords2d[, 1], fac)
      ylx <- split(object@stats$coords2d[, 2], fac)
      for(level in 1:nlev) {
        xbase <- object@stats$mean2d[level, 1]
        ybase <- object@stats$mean2d[level, 2]
        xlev <- xlx[[level]]
        ylev <- ylx[[level]]
        panel.segments(
          x0 = xbase, y0 = ybase,
          x1 =  xbase + object@g.args$starSize * (xlev - xbase),
          y1 =  ybase + object@g.args$starSize * (ylev - ybase),
          lty = plines$lty[level], lwd = plines$lwd[level], col = plines$col[level])
      }
    }
    
    ## draw points
    npoints <- nrow(object@stats$coords2d)
    ppoints <- object@adeg.par$ppoints
    if(length(fac) > 1) {
      ppoints <- lapply(object@adeg.par$ppoints, function(x, fac) {
        if(length(x) > length(fac))
          return(x)
        else {
          xlev <-  rep(x, length.out = nlev)
          return(xlev[fac])
        }
      }, fac = fac)
    }
    
    panel.points(x = object@stats$coords2d[, 1], y = object@stats$coords2d[, 2], type = "p", pch = ppoints$pch, cex = ppoints$cex, col = ppoints$col, alpha = ppoints$alpha, fill = ppoints$fill)
    
    ## draw labels
    if(any(object@adeg.par$plabels$cex > 0)) {
      center <- object@stats$mean2d
      adeg.panel.label(x = center[, 1], y = center[, 2] , labels = labels,  object@adeg.par$plabels)
    }
  })


triangle.class <- function(dfxyz, fac, wt = rep(1, NROW(fac)), labels = levels(fac), col = NULL, ellipseSize = 1, starSize = 1, chullSize = NULL, adjust = TRUE, 
  min3d = NULL, max3d = NULL, showposition = TRUE, facets = NULL, plot = TRUE, storeData = TRUE, add = FALSE, pos = -1, ...) {
  ## dfxyz: matrix/data.frame with 3 columns
  ## min3d, max3d: limits by default: c(0,0,0), c(1,1,1)
  
  thecall <- .expand.call(match.call())
  
  ## parameters sorted
  sortparameters <- sortparamADEg(...)

  ## facets
  if(!is.null(facets)) {
    if(NCOL(fac) == 1)
      object <- multi.facets.Tr(thecall, samelimits = sortparameters$g.args$samelimits)
    else 
      stop("Facets are not allowed multiple fac")
  }
  
  ## multiple fac
  else if(NCOL(fac) > 1) {
    object <- multi.variables.Tr(thecall, "fac")
  }
  
  ## simple ADEg graphic
  else {
    if(length(sortparameters$rest))
      warning(c("Unused parameters: ", paste(unique(names(sortparameters$rest)), " ", sep = "")), call. = FALSE)
    
     ## creation of the ADEg object
    g.args <- c(sortparameters$g.args, list(adjust = adjust, min3d = min3d, max3d = max3d, ellipseSize = ellipseSize, starSize = starSize, chullSize = chullSize, col = col))
    if(storeData)
    	tmp_data <- list(dfxyz = dfxyz, fac = fac, wt = wt, labels = labels, frame = sys.nframe() + pos, storeData = storeData)
    else
      tmp_data <- list(dfxyz = thecall$dfxyz, fac = thecall$fac, wt = thecall$wt, labels = thecall$labels, frame = sys.nframe() + pos, storeData = storeData)
    object <- new(Class = "Tr.class", data  = tmp_data, adeg.par = sortparameters$adepar, trellis.par = sortparameters$trellis, g.args = g.args, Call = match.call())
    
    ## preparation
    prepare(object)
    setlatticecall(object)
    if(showposition & add) {
      print("cannot show position and add") ## can be done, but modifies the meaning of the superposition
      showposition <- FALSE 
    }
    if(showposition)
      object <- new(Class = "ADEgS", ADEglist = list("triangle" = object, "positions" = .showpos(object)), positions = rbind(c(0, 0, 1, 1), c(0, 0.7, 0.3, 1)), add = matrix(0, ncol = 2, nrow = 2), Call = match.call())
    if(add)
      object <- add.ADEg(object)
  }

  if(!add & plot)
    print(object)
  invisible(object)
}

