#########################################################
## C1.gauss: here assumption: gaussian distribution   ###
#########################################################

setClass(
  Class = "C1.gauss",
  contains = "ADEg.C1"
)


setMethod(
  f = "initialize",
  signature = "C1.gauss",
  definition = function(.Object, data = list(score = NULL, fac = NULL, wt = NULL, frame = 0, storeData = TRUE), ...) {
    .Object <- callNextMethod(.Object, data = data, ...) ## ADEg.C1 initialize
    .Object@data$fac <- data$fac
    .Object@data$wt <- data$wt
    return(.Object)
  })


setMethod(
  f = "prepare",
  signature = "C1.gauss",
  definition = function(object) {
    nameobj <- deparse(substitute(object))
    
    ## pre-management of graphics parameters      
    oldparamadeg <- adegpar()
    on.exit(adegpar(oldparamadeg))
    adegtot <- adegpar(object@adeg.par)
    
    if(object@data$storeData) {
      fac <- object@data$fac
      score <- object@data$score
      wt <- object@data$wt
    } else {
      fac <- eval(object@data$fac, envir = sys.frame(object@data$frame))
      score <- eval(object@data$score, envir = sys.frame(object@data$frame))
      wt <- eval(object@data$wt, envir = sys.frame(object@data$frame))
    }
    nlev <- nlevels(as.factor(fac))
     
    ## If axes are plotted, put a label for axis
    if(adegtot$paxes$draw) {
      if(is.null(object@g.args$xlab) & !adegtot$p1d$horizontal)
        object@g.args$xlab <- "density"
      if(is.null(object@g.args$ylab) & adegtot$p1d$horizontal)
        object@g.args$ylab <- "density"
    }
    
    if(is.logical(object@g.args$col)) {
      if(object@g.args$col)
        adegtot$plabels$col <- adegtot$plabels$boxes$col <- adegtot$plines$col <- adegtot$ppolygons$col <- adegtot$ppolygons$border <- adegtot$ppalette$quali(nlev) 
    } else
      adegtot$plabels$col <- adegtot$plabels$boxes$col <- adegtot$plines$col <- adegtot$ppolygons$col <- adegtot$ppolygons$border <- rep(object@g.args$col, length.out = nlev)
    
    ## if fill is FALSE, polygons density curves are transparent
    if(!object@g.args$fill)
      adegtot$ppolygons$col <- "transparent"
    
    ## object modification before calling inherited method
    object@adeg.par <- adegtot
    callNextMethod() ## prepare graph
   
    ## statistics calculus
    object@stats$means <- meanfacwt(score, fac, wt)
    object@stats$var <- varfacwt(score, fac)
    
    ## here steps fixed, could be a argument of s1d.gauss
    steps <- object@g.args$steps
    nind <- table(fac)
    gausscurv <- list()
    if(object@adeg.par$p1d$horizontal)
      xx <- seq(from = object@g.args$xlim[1], to = object@g.args$xlim[2], length.out = steps)
    else
      xx <- seq(from = object@g.args$ylim[1], to = object@g.args$ylim[2], length.out = steps)
    
    for(i in 1:nlev) {
      if(nind[i] == 0)
        gausscurv[[i]] <- NA
      else
        gausscurv[[i]] <- dnorm(xx, mean = object@stats$means[i], sd = sqrt(object@stats$var[i]))
    }
    names(gausscurv) <- levels(fac)
    
    lead <- ifelse(object@adeg.par$p1d$reverse, 1 , -1)
    
    if(object@adeg.par$p1d$horizontal && is.null(object@g.args$ylim))
      object@g.args$ylim <- c(0, max(sapply(gausscurv, FUN = function(x) {ifelse(is.na(x[1]), 0, max(x)) / 0.85})))
    if(object@adeg.par$p1d$horizontal) {
    	ref <- ifelse(object@adeg.par$p1d$reverse, 2, 1)
      margin <- object@g.args$ylim[ref]
      if(object@adeg.par$p1d$rug$draw)
        margin <- object@adeg.par$p1d$rug$margin * abs(diff(object@g.args$ylim))
      object@s.misc$rug <- object@g.args$ylim[ref]
      object@g.args$ylim[ref] <- object@g.args$ylim[ref] + lead * margin
    }
    
    if(!object@adeg.par$p1d$horizontal && is.null(object@g.args$xlim))
      object@g.args$xlim <- c(0, max(sapply(gausscurv, FUN = function(x) {ifelse(is.na(x[1]), 0, max(x)) / 0.85})))
    if(!object@adeg.par$p1d$horizontal) {
      ref <- ifelse(object@adeg.par$p1d$reverse, 2, 1)
      margin <- object@g.args$xlim[ref]
      if(object@adeg.par$p1d$rug$draw)
        margin <- object@adeg.par$p1d$rug$margin * abs(diff(object@g.args$xlim))
      object@s.misc$rug <- object@g.args$xlim[ref]
      object@g.args$xlim[ref] <- object@g.args$xlim[ref] + lead * margin
    }
    
    object@stats$gausscurves <- gausscurv
    assign(nameobj, object, envir = parent.frame())
  })


setMethod(
  f = "panel",
  signature = "C1.gauss",
  definition = function(object, x, y) {
    ## Drawing gauss curves as polygons (filled or not)
    ## one polygon per level
    ## y is the score

    ## get some parameters
    pscore <- object@adeg.par$p1d
    curvess <- object@stats$gausscurves
    labels <- names(curvess)
    lims <- current.panel.limits(unit = "native")
    
    if(object@data$storeData)
      fac <- object@data$fac
    else
      fac <- eval(object@data$fac, envir = sys.frame(object@data$frame))
    nlev <- nlevels(as.factor(fac))
       
    ppoly <- lapply(object@adeg.par$ppolygons, FUN = function(x) rep(x, length.out = nlev))
    plabels <- lapply(object@adeg.par$plabels, FUN = function(x) rep(x, length.out = nlev))

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
    
    ## Starts the display
    ## depends on the parametres horizontal and reverse
    lead <- ifelse(pscore$reverse, -1, 1)   
    if(pscore$horizontal) {
      ## horizontal drawing
      margin <- 0
      xx <- seq(from = lims$xlim[1], to = lims$xlim[2], length.out = object@g.args$steps)
      if(pscore$rug$draw)
        margin <- if(is.unit(pscore$rug$margin)) convertUnit(pscore$rug$margin, typeFrom = "dimension", unitTo = "native", axisFrom = "y", valueOnly = TRUE) else pscore$rug$margin
      margin <- ifelse(pscore$reverse, lims$ylim[2], lims$ylim[1]) + lead * margin
       
      for(i in 1:nlev) {
        if(!is.na(curvess[[i]][1])) {
          y <- margin + lead * curvess[[i]]
          panel.polygon(x = c(lims$xlim[1], xx, lims$xlim[2]), y = c(margin, y, margin) , border = ppoly$border[i],
                        col = ppoly$col[i], lty = ppoly$lty[i], lwd = ppoly$lwd[i], alpha = ppoly$alpha[i])
          if(nlev > 1) {
            ## indicate levels names for each curve
            ymaxindex <- which.max(curvess[[i]]) ## places at the maximum
            panel.text(x = xx[ymaxindex], y = y[ymaxindex], labels = names(curvess)[i], pos = ifelse(pscore$reverse, 1, 3), col = plabels$col[i],
                       cex = plabels$cex[i], alpha = plabels$alpha[i], srt = srt)
          }
        }
      }
    } else {
      ## vertical drawing
      margin <- 0
      yy <- seq(from = lims$ylim[1], to = lims$ylim[2], length.out = object@g.args$steps)
      if(pscore$rug$draw)
        margin <- if(is.unit(pscore$rug$margin)) convertUnit(pscore$rug$margin, typeFrom = "dimension", unitTo = "native", axisFrom = "x", valueOnly = TRUE) else pscore$rug$margin
      margin <- ifelse(pscore$reverse, lims$xlim[2], lims$xlim[1]) + lead * margin
      for(i in 1:nlev) {
        if(!is.na(curvess[[i]][1])) {
          x <- margin + lead * curvess[[i]]
          panel.polygon(x = c(margin, x, margin), y = c(lims$ylim[1], yy, lims$ylim[2]), border = ppoly$border[i],
                        col = ppoly$col[i], lty = ppoly$lty[i], lwd = ppoly$lwd[i], alpha = ppoly$alpha[i])
          if(nlev > 1) {
            xmaxindex <- which.max(curvess[[i]])
            panel.text(x = x[xmaxindex], y = yy[xmaxindex], labels = names(curvess)[i], col = plabels$col[i], pos = ifelse(pscore$reverse, 2, 4),
                       cex = plabels$cex[i], alpha = plabels$alpha[i], srt = srt)
          }
        }
      }
    }
  })


s1d.gauss <- function(score, fac = gl(1, NROW(score)), wt = rep(1, NROW(score)), steps = 200, col = TRUE, fill = TRUE, facets = NULL,
                      plot = TRUE, storeData = TRUE, add = FALSE, pos = -1, ...) {
                      
  thecall <- .expand.call(match.call())
  
  ## parameters management
  sortparameters <- sortparamADEg(...)
  
  ## facets
  if(!is.null(facets)) {
    if(NCOL(score) == 1 & NCOL(fac) == 1)
      object <- multi.facets.C1(thecall, sortparameters$adepar, samelimits = sortparameters$g.args$samelimits)
    else 
      stop("Facets are not allowed with multiple scores and/or multiple fac")
  }
  
  ## multiple scores
  else if(NCOL(score) > 1) {
    if(NCOL(fac) == 1)
      object <- multi.score.C1(thecall)
    else 
      stop("Multiple scores are not allowed with multiple fac")
  }
  
  ## multiple fac
  else if(NCOL(fac) > 1) {
    object <- multi.variables.C1(thecall, "fac")
  }
  
  ## simple ADEg graphic
  else {
    if(length(sortparameters$rest))
      warning(c("Unused parameters: ", paste(unique(names(sortparameters$rest)), " ", sep = "")), call. = FALSE)
    
    ## creation of the ADEg object
    g.args <- c(sortparameters$g.args, list(steps = steps, fill = fill, col = col))
    if(storeData)
    	tmp_data <- list(score = score, fac = fac, wt = wt, frame = sys.nframe() + pos, storeData = storeData)
    else
      tmp_data <- list(score = thecall$score, fac = thecall$fac, wt = thecall$wt, frame = sys.nframe() + pos, storeData = storeData)
    object <- new(Class = "C1.gauss", data = tmp_data, adeg.par = sortparameters$adepar, trellis.par = sortparameters$trellis, g.args = g.args, Call = match.call())

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
