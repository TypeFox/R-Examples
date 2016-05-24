##########################################################################
##                            s.corcircle                               ##
##########################################################################

setClass(
  Class = "S2.corcircle",
  contains = "ADEg.S2",
)


setMethod(
  f = "initialize",
  signature = "S2.corcircle",
  definition = function(.Object, data = list(dfxy = NULL, xax = 1, yax = 2, labels = NULL, frame = 0, storeData = TRUE), ...) {
    .Object <- callNextMethod(.Object, data = data, ...) ## ADEg.S2 initialize
    .Object@data$labels <- data$labels
    return(.Object)
  })


setMethod(
  f = "prepare",
  signature = "S2.corcircle",
  definition = function(object) {
    name_obj <- deparse(substitute(object))
    
    ## pre-management of graphics parameters      
    oldparamadeg <- adegpar()
    on.exit(adegpar(oldparamadeg))
    adegtot <- adegpar(object@adeg.par)
    
    ## prepare grid
    getgrid <- function(nbgrid = 10) {
      cgrid <- signif(2 / nbgrid, 2)
      h0 <- c(rev(seq(0, -1, by = -cgrid)), seq(0 + cgrid, 1, by = cgrid)) ## force that 0 is represented by the grid
      cgrid <- diff(h0)[1]
      coord <- rep(0, length(h0))
      for(i in 1:length(h0))
        coord[i] <- sqrt(1 - h0[i] * h0[i])
      return(list(x0 = c(h0, -coord), x1 = c(h0, coord), y0 = c(-coord, h0), y1 = c(coord, h0), d = cgrid))
    }
    if(adegtot$pgrid$draw || adegtot$paxes$draw)
      object@s.misc$backgrid <- getgrid(adegtot$pgrid$nint)      
    if(is.null(object@adeg.par$ppoints$cex))
      adegtot$ppoints$cex <- 0
    
    ## object modification before calling inherited method
    object@adeg.par <- adegtot
    callNextMethod() ## prepare graph
    assign(name_obj, object, envir = parent.frame())
  })


setMethod(
  f = "panel",
  signature = "S2.corcircle",
  definition = function(object, x, y) {
    panel.arrows(x0 = 0, y0 = 0, y1 = y, x1 = x, angle = object@adeg.par$parrows$angle,
                 length = object@adeg.par$parrows$length, ends = object@adeg.par$parrows$end,
                 lwd = object@adeg.par$plines$lwd, col = object@adeg.par$plines$col, lty = object@adeg.par$plines$lty)
    ## labels and boxes           
    plabels <- object@adeg.par$plabels
    
    pos <- .textpos(x, y, origin = c(0, 0))
    
    if(object@data$storeData)
      labels <- object@data$labels
    else
      labels <- eval(object@data$labels, envir = sys.frame(object@data$frame))
       
    test <- .textsize(labels, plabels)
    w <- test$w
    h <- test$h
    plabels$optim <- FALSE ## no optimisation
    adeg.panel.label(x = x + pos[1, ] * w / 2, y = y + pos[2, ] * h / 2, labels = labels, plabels = plabels)
  })


s.corcircle <- function(dfxy, xax = 1, yax = 2, labels = row.names(as.data.frame(dfxy)), fullcircle = TRUE,
                      	facets = NULL, plot = TRUE, storeData = TRUE, add = FALSE, pos = -1, ...) {
                        
  ## evaluation of some parameters (required for multiplot)
  thecall <- .expand.call(match.call())
  df <- try(as.data.frame(eval(thecall$dfxy, envir = sys.frame(sys.nframe() + pos))), silent = TRUE)
  if((class(df) == "try-error") | is.null(thecall$dfxy)) ## non convenient dfxy argument
    stop("non convenient selection for dfxy (can not be converted to dataframe)")
  
  ## parameters sorted
  sortparameters <- sortparamADEg(...)
  
  ## facets
  if(!is.null(facets)) { 
    if((length(xax) == 1 & length(yax) == 1))
      object <- multi.facets.S2(thecall, sortparameters$adepar)
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
    g.args <- c(sortparameters$g.args, list(fullcircle = fullcircle))
    if(storeData)
      tmp_data <- list(dfxy = dfxy, xax = xax, yax = yax, labels = labels, frame = sys.nframe() + pos, storeData = storeData)
    else
    	tmp_data <- list(dfxy = thecall$dfxy, xax = xax, yax = yax, labels = thecall$labels, frame = sys.nframe() + pos, storeData = storeData)
    object <- new(Class = "S2.corcircle", data = tmp_data, adeg.par = sortparameters$adepar, trellis.par = sortparameters$trellis, g.args = g.args, Call = as.call(thecall))
    
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
