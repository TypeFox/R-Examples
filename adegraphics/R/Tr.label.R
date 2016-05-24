###############################################
##             triangle.label                ##
###############################################

setClass(
  Class = "Tr.label",
  contains = "ADEg.Tr"
)


setMethod(
  f = "initialize",
  signature  = "Tr.label",
  definition = function(.Object, data = list(dfxyz = NULL, labels = NULL, frame = 0, storeData = TRUE), ...) {
    .Object <- callNextMethod(.Object, data = data, ...) ## ADEg.Tr initialize
    .Object@data$labels <- data$labels
    return(.Object)
  })


setMethod(
  f = "prepare",
  signature = "Tr.label",
  definition = function(object) {
    name_obj <- deparse(substitute(object))
    
    if(object@data$storeData) {
      labels <- object@data$labels
      df <- object@data$dfxyz
    } else {
      labels <- eval(object@data$labels, envir = sys.frame(object@data$frame))
      df <- eval(object@data$dfxyz, envir = sys.frame(object@data$frame))
    }
    
    ## pre-management of graphics parameters      
    oldparamadeg <- adegpar()
    on.exit(adegpar(oldparamadeg))
    adegtot <- adegpar(object@adeg.par)
    
    ## change some parameter values
    if((is.null(object@adeg.par$plabels$boxes$draw) & adegtot$plabels$optim) | (is.null(object@adeg.par$plabels$boxes$draw) & length(labels) > 1000))  
      adegtot$plboxes$draw <- FALSE
    
    if(object@g.args$addmean) {
      default <- list(pch = 20, col = "black", cex = 2)
      if(is.list(object@g.args$meanpar))
        object@g.args$meanpar <- modifyList(default, object@g.args$meanpar, keep.null = TRUE)
      else {
        if(!is.null(object@g.args$meanpar))
          stop("meanpar must be a list of graphical parameters (pch, col, cex)", call. = FALSE)
        else
          object@g.args$meanpar <- default
      }
    }
    
    if(object@g.args$addaxes | object@g.args$addmean) {
      ## lines (axes or mean)
      default <- list(col = "black", lwd = 1, lty = 1)
      if(is.list(object@g.args$axespar))
        object@g.args$axespar <- modifyList(default, object@g.args$axespar, keep.null = TRUE)
      else {
        if(!is.null(object@g.args$axespar))
          stop("axespar must be a list of graphical parameters (lwd, col, lty)", call. = FALSE)
        else
          object@g.args$axespar <- default
      }
      
      ## point (axes or mean)
      default <- list(pch = 20, col = "black", cex = 2)
      if(is.list(object@g.args$meanpar))
        object@g.args$meanpar <- modifyList(default, object@g.args$meanpar, keep.null = TRUE)
      else {
        if(!is.null(object@g.args$meanpar))
          stop("meanpar must be a list of graphical parameters (pch, col, cex)", call. = FALSE)
        else
          object@g.args$meanpar <- default
      }
    }
    
    ## object modification before calling inherited method
    object@adeg.par <- adegtot
    callNextMethod() ## prepare graph
    
    ## calculate 2D coordinates
    df <- sweep(df, 1, rowSums(df), "/")
    object@stats$coords2d <- .coordtotriangleM(df, mini3 = object@g.args$min3d, maxi3 = object@g.args$max3d)[, 2:3]
    
    assign(name_obj, object, envir = parent.frame())
  })


setMethod(
  f = "panel",
  signature = "Tr.label",
  definition = function(object, x, y) {            
    
    if(object@data$storeData) {
      labels <- object@data$labels
      df <- object@data$dfxyz
    } else {
      labels <- eval(object@data$labels, envir = sys.frame(object@data$frame))
      df <- eval(object@data$dfxyz, envir = sys.frame(object@data$frame))
    }
    
    ## draw points and labels
    if(any(object@adeg.par$ppoints$cex > 0))
      panel.points(object@stats$coords2d[, 1], object@stats$coords2d[, 2], pch = object@adeg.par$ppoints$pch, cex = object@adeg.par$ppoints$cex, col = object@adeg.par$ppoints$col, alpha = object@adeg.par$ppoints$alpha, fill = object@adeg.par$ppoints$fill)
    if(any(object@adeg.par$plabels$cex > 0))
      adeg.panel.label(object@stats$coords2d[, 1], object@stats$coords2d[, 2], labels, object@adeg.par$plabels)
    
    ## addmean or addaxes
    if(object@g.args$addmean | object@g.args$addaxes) {
      df <- sweep(df, 1, rowSums(df), "/")
      mini3 <- object@g.args$min3d
      maxi3 <- object@g.args$max3d
      
      m3 <- colMeans(df)
      mxy <- .coordtotriangleM(t(as.matrix(m3)), mini3 = mini3, maxi3 = maxi3)[-1]
      if(object@g.args$addmean) {
        ## axis points: putting means on the axis A
        axp3 <- rbind(c(m3[1], mini3[2], 1 - m3[1] - mini3[2]),
                      c(1 - m3[2] -mini3[3], m3[2], mini3[3]),
                      c(mini3[1], 1 - m3[3] - mini3[1], m3[3]))
        axpxyz <- .coordtotriangleM(axp3, mini3 = mini3, maxi3 = maxi3)
        
        ## drawing lines for means
        apply(axpxyz, 1, FUN = function(x) {
          do.call("panel.lines", c(list(x = c(x[2], mxy[1]), y = c(x[3], mxy[2])), object@g.args$axespar))
        })
        do.call("panel.points", c(list(x = c(mxy[1], axpxyz[, 2]), y = c(mxy[2], axpxyz[, 3])), object@g.args$meanpar))
        panel.text(x = axpxyz[, 2], y = axpxyz[, 3], labels = as.character(round(m3, digits = 4)), pos = c(2, 1, 4))
      }
      
      if(object@g.args$addaxes) {
        axx <- dudi.pca(df, scale = FALSE, scannf = FALSE)$c1
        cornerp <- object@s.misc$cornerp
        a1 <- axx[, 1]
        x1 <- a1[1] * cornerp$A + a1[2] * cornerp$B + a1[3] * cornerp$C
        do.call("panel.segments", c(list(x0 = mxy[1] - x1[1], x1 = mxy[1] +  x1[1], y0 = mxy[2] - x1[2], y1 = mxy[2] + x1[2]), object@g.args$axespar))
        a2 <- axx[, 2]
        x1 <- a2[1] * cornerp$A + a2[2] * cornerp$B + a2[3] * cornerp$C
        do.call("panel.segments", c(list(x0 = mxy[1] - x1[1], x1 = mxy[1] + x1[1], y0 = mxy[2] - x1[2], y1 = mxy[2] + x1[2]), object@g.args$axespar))
        do.call("panel.points", c(list(x = mxy[1], y = mxy[2]), object@g.args$meanpar))
      }
    }
  })


triangle.label <- function(dfxyz, labels = rownames(dfxyz), adjust = TRUE, min3d = NULL, max3d = NULL, addaxes = FALSE, addmean = FALSE, meanpar = NULL, axespar = NULL, 
  												 showposition = TRUE, facets = NULL, plot = TRUE, storeData = TRUE, add = FALSE, pos = -1, ...) {
  ## dfxyz: matrix/data.frame with 3 columns
  ## min3d, max3d: limits by default: c(0,0,0), c(1,1,1)
  ## addaxes: should we draw pca axes
  ## addmean: should we draw mean
  
  thecall <- .expand.call(match.call())
  
  ## parameters sorted
  sortparameters <- sortparamADEg(...)
  
  ## facets
  if(!is.null(facets)) {
    object <- multi.facets.Tr(thecall, samelimits = sortparameters$g.args$samelimits)
  }
  
  ## simple ADEg graphic
  else {
    if(length(sortparameters$rest))
      warning(c("Unused parameters: ", paste(unique(names(sortparameters$rest)), " ", sep = "")), call. = FALSE)
    
    ## creation of the ADEg object
    g.args <- c(sortparameters$g.args, list(adjust = adjust, min3d = min3d, max3d = max3d, addaxes = addaxes, addmean = addmean, meanpar = meanpar, axespar = axespar))
    if(storeData)
    	tmp_data <- list(dfxyz = dfxyz, labels = labels, frame = sys.nframe() + pos, storeData = storeData)
    else
      tmp_data <- list(dfxyz = thecall$dfxyz, labels = thecall$labels, frame = sys.nframe() + pos, storeData = storeData)
    object <- new(Class = "Tr.label", data = tmp_data, adeg.par = sortparameters$adepar, trellis.par = sortparameters$trellis, g.args = g.args, Call = match.call())
    
    ## preparation
    prepare(object)
    setlatticecall(object)
    if(showposition & add) {
      warning("cannot show position and add") ## can be done, but modifies the meaning of the superposition
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

