##########################################################################
##                            s.density                                 ##
##########################################################################

setClass(
  Class = "S2.density",
  contains = "ADEg.S2"
)


## no initialize function (use ADEg.S2 by inheritance)
setMethod(
  f = "prepare",
  signature = "S2.density",
  definition = function(object) {
    name_obj <- deparse(substitute(object))
    
    ## pre-management of graphics parameters      
    oldparamadeg <- adegpar()
    on.exit(adegpar(oldparamadeg))
    adegtot <- adegpar(object@adeg.par)
    
    if(object@data$storeData)
      dfxy <- object@data$dfxy
    else 
      dfxy <- eval(object@data$dfxy, envir = sys.frame(object@data$frame))
    
    if(is.null(object@adeg.par$plabels))
      adegtot$plabels$cex <- 0
    if(is.null(object@adeg.par$porigin$include) & (any(names(object@g.args) %in% c("Sp", "nbobject"))))
      adegtot$porigin$include <- FALSE
     
    ## object modification before calling inherited method
    object@adeg.par <- adegtot
    callNextMethod() ## prepare graph
    
    ## compute density using bkde2D (KernSmooth package)    
    ## bandwidth and gridsize can be provided by user. range.x allows computation for all the panel (even with no points)
    if(is.null(object@g.args$bandwidth))
      object@g.args$bandwidth <- diff(apply(dfxy[, c(object@data$xax, object@data$yax)], 2, quantile, probs = c(0.05, 0.95), na.rm = TRUE)) / 25
    if(min(object@g.args$bandwidth) <= 0)
      stop("'bandwidth' must be strictly positive")
    object@g.args$threshold <- min(max(0, object@g.args$threshold), 1)
    
    object@stats$densit <- bkde2D(dfxy[, c(object@data$xax[1], object@data$yax[1])], bandwidth = object@g.args$bandwidth, gridsize = rep(object@g.args$gridsize, length.out = 2))
    ## TODO: as in s.image, remove points (only) where density is null
    ## use expand.grid....           
    assign(name_obj, object, envir = parent.frame())
  })


setMethod(
    f = "panel",
    signature = "S2.density",
    definition = function(object, x, y) {
        densit <- object@stats$densit
        if(is.null(object@g.args$col))
            col <- object@adeg.par$ppalette$quanti(255)
        else
            col <- object@g.args$col
        
        transformation <- function(x) x
        densityy <- array(transformation(densit$fhat), dim = dim(densit$fhat))            
        if(object@g.args$region)
            panel.levelplot(x = rep(densit$x1, length(densit$x2)),
                            y = rep(densit$x2, each = length(densit$x1)),
                            z = densityy,
                            at = c(-.Machine$double.eps, seq(from = max(densit$fhat) * object@g.args$threshold + .Machine$double.eps,
                                to = 1.01 * max(densit$fhat), length = length(col) + 2)),
                            col.regions = c("transparent", col),
                            subscripts = TRUE)
        if(object@g.args$contour)
            panel.levelplot(x = rep(densit$x1, length(densit$x2)), 
                            y = rep(densit$x2, each = length(densit$x1)), 
                            z = densityy,
                            labels = object@adeg.par$plabels, 
                            label.style = if(object@adeg.par$plabels$srt == "horizontal") "flat" else "align",  ## also exist "mixed" not used here
                            at = c(-.Machine$double.eps, seq(from = max(densit$fhat) * object@g.args$threshold + .Machine$double.eps, 
                                to = 1.01 * max(densit$fhat), length = object@g.args$nclass + 1)), 
                            col.regions = c("transparent", col), 
                            subscripts = TRUE, 
                            region = FALSE, 
                            contour = TRUE)
        
        
        ## show nrpoints outilers
        if(object@g.args$nrpoints > 0) {
            ## copy of panel.smoothScatter
            ixm <- round((x - densit$x1[1]) / (densit$x1[length(densit$x1)] - densit$x1[1]) * (length(densit$x1) - 1))
            iym <- round((y - densit$x2[1]) / (densit$x2[length(densit$x2)] - densit$x2[1]) * (length(densit$x2) - 1))
            idens <- densityy[1 + iym * length(densit$x1) + ixm]
            nrpoints <- min(nrow(x), ceiling(object@g.args$nrpoints))
            sel <- order(idens, decreasing = FALSE)[1:nrpoints]
            
            panel.points(x[sel], y[sel], pch = object@adeg.par$ppoints$pch, cex = object@adeg.par$ppoints$cex, col = object@adeg.par$ppoints$col, fill = object@adeg.par$ppoints$fill)
        }
    })


s.density <- function(dfxy, xax = 1, yax = 2, bandwidth = NULL, gridsize = c(450L, 450L), nrpoints = 300, threshold = 0.1, col = NULL, 
  contour = FALSE, region = !contour, nclass = 8, facets = NULL, plot = TRUE, storeData = TRUE, add = FALSE, pos = -1, ...) {
  
  ## evaluation of some parameters
  thecall <- .expand.call(match.call())
  df <- try(as.data.frame(eval(thecall$dfxy, envir = sys.frame(sys.nframe() + pos))), silent = TRUE)
  if((class(df) == "try-error") | is.null(thecall$dfxy)) ## non convenient dfxy argument
    stop("non convenient selection for dfxy (can not be converted to dataframe)")

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
    g.args <- c(sortparameters$g.args, list(bandwidth = bandwidth, gridsize = gridsize, threshold = threshold, col = col, nrpoints = nrpoints, contour = contour, region = region, nclass = nclass))
    if(storeData)
    	tmp_data <- list(dfxy = dfxy, xax = xax, yax = yax, frame = sys.nframe() + pos, storeData = storeData)
    else
      tmp_data <- list(dfxy = thecall$dfxy, xax = xax, yax = yax, frame = sys.nframe() + pos, storeData = storeData)
    object <- new(Class = "S2.density", data = tmp_data, adeg.par = sortparameters$adepar, trellis.par = sortparameters$trellis, g.args = g.args, Call = as.call(thecall))
    
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
