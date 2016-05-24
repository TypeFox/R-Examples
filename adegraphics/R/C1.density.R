#####################################################################
## S1.density to compare with S1.gauss afterwards                  ##
## TODO: reverse/vertical mettre a l'echelle distribution calculee ##
## Dans l'idée S1.density plutot si pas de factor...               ##
#####################################################################

setClass(
  Class = "C1.density",
  contains = "ADEg.C1"
)


setMethod(
  f = "initialize",
  signature = "C1.density",
  definition = function(.Object, data = list(score = NULL, fac = NULL, frame = 0, storeData = TRUE), ...) {
    .Object <- callNextMethod(.Object, data = data, ...) ## ADEg.C1 initialize
    .Object@data$fac <- data$fac
    return(.Object)
  })


### densities calculations according to user parameters and score/factor
setMethod(
  f = "prepare",
  signature = "C1.density",
  definition = function(object) {
    nameobj <- deparse(substitute(object))
    
    ## pre-management of graphics parameters      
    oldparamadeg <- adegpar()
    on.exit(adegpar(oldparamadeg))
    adegtot <- adegpar(object@adeg.par)
    
    if(object@data$storeData) {
      fac <- object@data$fac
      score <- object@data$score
    } else {
      fac <- eval(object@data$fac, envir = sys.frame(object@data$frame))
      score <- eval(object@data$score, envir = sys.frame(object@data$frame))
    }
    score <- as.matrix(score)[, 1]  ## to manage 'score' when it is a data.frame with only one column
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
    
    scores <- split(score, fac)
    densit <- vector(mode = "list", length = length(scores))
    names(densit) <- names(scores)
    ## estimate density for each level of the factor
    for(i in 1:length(scores)) {
      if(length(scores[[i]]) == 0) {
        ## no data in the given level
        densit[[i]] <- list(x = NA, y = NA)
      } else {
        if(!is.null(object@g.args$bandwidth))
          densit[[i]] <- bkde(scores[[i]], kernel = object@g.args$kernel, bandwidth = object@g.args$bandwidth, gridsize = object@g.args$gridsize)
        else
          densit[[i]] <- bkde(scores[[i]], kernel = object@g.args$kernel, gridsize = object@g.args$gridsize)        
      }
    }
   
    lead <- ifelse(object@adeg.par$p1d$reverse, 1 , -1)
    
    if(object@adeg.par$p1d$horizontal && is.null(object@g.args$ylim))
      object@g.args$ylim <- c(0, max(sapply(densit, FUN = function(x) {ifelse(is.na(x$y[1]), 0, max(x$y))}) / 0.85))
    if(object@adeg.par$p1d$horizontal) {
    	ref <- ifelse(object@adeg.par$p1d$reverse, 2, 1)
      margin <- object@g.args$ylim[ref]
      if(object@adeg.par$p1d$rug$draw)
        margin <- object@adeg.par$p1d$rug$margin * abs(diff(object@g.args$ylim))
      object@s.misc$rug <- object@g.args$ylim[ref]
      object@g.args$ylim[ref] <- object@g.args$ylim[ref] + lead * margin
    }
    
    if(!object@adeg.par$p1d$horizontal && is.null(object@g.args$xlim))
      object@g.args$xlim <- c(0, max(sapply(densit, FUN = function(x) {ifelse(is.na(x$y[1]), 0, max(x$y))}) / 0.85))
    if(!object@adeg.par$p1d$horizontal) {
      ref <- ifelse(object@adeg.par$p1d$reverse, 2, 1)
      margin <- object@g.args$xlim[ref]
      if(object@adeg.par$p1d$rug$draw)
        margin <- object@adeg.par$p1d$rug$margin * abs(diff(object@g.args$xlim))
      object@s.misc$rug <- object@g.args$xlim[ref]
      object@g.args$xlim[ref] <- object@g.args$xlim[ref] + lead * margin
    }
    
    object@stats$densit <- densit
    assign(nameobj, object, envir = parent.frame())
  })


setMethod(
  f = "panel",
  signature = "C1.density",
  definition = function(object, x, y) {
    ## Drawing densities as polygons (filled or not)
    ## one polygon per level
    ## y is the score

    ## get some parameters
    pscore <- object@adeg.par$p1d
    curvess <- object@stats$densit
    labels <- names(curvess)
    lims <- current.panel.limits(unit = "native")
   
    if(object@data$storeData)
      fac <- object@data$fac
    else
      fac <- eval(object@data$fac, envir = sys.frame(object@data$frame))
    nlev <- nlevels(as.factor(fac))
   
    ppoly <- lapply(object@adeg.par$ppolygons, FUN = function(x) rep(x, length.out = nlev))
    plabels <- lapply(object@adeg.par$plabels, FUN = function(x) rep(x, length.out = nlev))
    
    y <- split(y, fac)
        
    ## manage string rotation
    srt <- 0
    if(is.numeric(plabels$srt[1]))
      srt <- plabels$srt[1]
    else{
      if(plabels$srt[1] == "horizontal")
        srt <- 0
      else if(plabels$srt[1] == "vertical")
        srt <- 90
    }
 
    ##  Starts the display
    ## depends on the parametres horizontal and reverse
    lead <- ifelse(pscore$reverse, -1, 1)
    if(pscore$horizontal) {
      ## horizontal drawing
      margin <- 0
      if(pscore$rug$draw)
        margin <- if(is.unit(pscore$rug$margin)) convertUnit(pscore$rug$margin, typeFrom = "dimension", unitTo = "native", axisFrom = "y", valueOnly = TRUE) else pscore$rug$margin
      margin <- ifelse(pscore$reverse, lims$ylim[2], lims$ylim[1]) + lead * margin
      
      for(i in 1:nlev) {
        if(!is.na(curvess[[i]]$y[1])) {
          y <- margin + lead * curvess[[i]]$y
          panel.polygon(x = c(min(curvess[[i]]$x), curvess[[i]]$x, max(curvess[[i]]$y)), y = c(margin, y, margin), border = ppoly$border[i],
                        col = ppoly$col[i], lty = ppoly$lty[i], lwd = ppoly$lwd[i], alpha = ppoly$alpha[i])
          if(nlev > 1) {
            ## indicate levels names for each curve
            ymaxindex <- which.max(curvess[[i]]$y) ## places at the maximum
            panel.text(x = curvess[[i]]$x[ymaxindex], y = y[ymaxindex], labels = names(curvess)[i], pos = ifelse(pscore$reverse, 1, 3), col = plabels$col[i],
                       cex = plabels$cex[i], alpha = plabels$alpha[i], srt = srt)
          }
        }
      }
    } else {
      ## vertical drawing
      margin <- 0
      if(pscore$rug$draw)
        margin <- if(is.unit(pscore$rug$margin)) convertUnit(pscore$rug$margin, typeFrom = "dimension", unitTo = "native", axisFrom = "x", valueOnly = TRUE) else pscore$rug$margin
      margin <- ifelse(pscore$reverse, lims$xlim[2], lims$xlim[1]) + lead * margin
   
      for(i in 1:nlev) {
        if(!is.na(curvess[[i]]$y[1])) {
          x <- margin + lead * curvess[[i]]$y
          panel.polygon(x = c(margin, x, margin), y = c(min(curvess[[i]]$x), curvess[[i]]$x, max(curvess[[i]]$x)), border = ppoly$border[i],
                        col = ppoly$col[i], lty = ppoly$lty[i], lwd = ppoly$lwd[i], alpha = ppoly$alpha[i])
          if(nlev > 1) {
            ## indicate levels names for each curve
            xmaxindex <- which.max(curvess[[i]]$y)
            panel.text(x = x[xmaxindex], y = curvess[[i]]$x[xmaxindex], labels = names(curvess)[i], pos = ifelse(pscore$reverse, 2, 4), col = plabels$col[i], cex = plabels$cex[i], alpha = plabels$alpha[i], srt = srt)
          }
        }
      }      
    }        
  })


## s1d.density: user function
## kernel, bandwidth and gridsize directly passed to the bkde function (for density calculation)
## if fill is FALSE, polygons density curves are transparent
s1d.density <- function(score, fac = gl(1, NROW(score)), kernel = c("normal", "box", "epanech", "biweight", "triweight"), bandwidth = NULL, gridsize = 450, col = TRUE, fill = TRUE, facets = NULL,
                        plot = TRUE, storeData = TRUE, add = FALSE, pos = -1, ...) {

	thecall <- .expand.call(match.call())
  
  ## parameters sorted
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
    g.args <- c(sortparameters$g.args, list(kernel = match.arg(kernel), bandwidth = bandwidth, gridsize = gridsize, fill = fill, col = col))
    if(storeData)
    	tmp_data <- list(score = score, fac = fac, frame = sys.nframe() + pos, storeData = storeData)
    else
      tmp_data <- list(score = thecall$score, fac = thecall$fac, frame = sys.nframe() + pos, storeData = storeData)
    object <- new(Class = "C1.density", data = tmp_data, adeg.par = sortparameters$adepar, trellis.par = sortparameters$trellis, g.args = g.args, Call = match.call())

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
