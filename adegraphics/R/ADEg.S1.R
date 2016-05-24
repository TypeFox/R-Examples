####################################################
##               Uni-dimensionnal plot            ##
####################################################

setClass(
  Class = "ADEg.S1",
  contains = c("ADEg", "VIRTUAL"),
  slots = c(data = "list")
  )

  
setMethod(
  f = "initialize",
  signature = "ADEg.S1",
  definition = function(.Object, data = list(score = NULL, at = NULL, frame = 0, storeData = TRUE), ...) {
    .Object <- callNextMethod(.Object, ...) ## ADEg initialize
    .Object@data <- data
    return(.Object)
  })


## prepare: grid calculations
## reset limits and sets axis information for lattice
setMethod(
  f = "prepare",
  signature = "ADEg.S1",
  definition = function(object) {
    name_obj <- deparse(substitute(object))
    if(object@data$storeData) {
      score <- object@data$score
      at <- object@data$at
    } else {
      score <- eval(object@data$score, envir = sys.frame(object@data$frame))
      at <- eval(object@data$at, envir = sys.frame(object@data$frame))
    }
    
    if(inherits(object, "S1.boxplot")){
      if(object@data$storeData) {
        fac <- object@data$fac
      } else {
        fac <- eval(object@data$fac, envir = sys.frame(object@data$frame))
      }
    }
    
    score <- as.matrix(score)[, 1]  ## to manage 'score' when it is a data.frame with only one column
    
    ## limits and scale
    minX <- min(score)
    maxX <- max(score)
    if(object@adeg.par$p1d$horizontal & !is.null(object@g.args$xlim) & is.null(object@s.misc$hori.update)) {
      minX <- object@g.args$xlim[1]
      maxX <- object@g.args$xlim[2]
    }

    if(!object@adeg.par$p1d$horizontal & !is.null(object@g.args$ylim) & is.null(object@s.misc$hori.update)) {
      minX <- object@g.args$ylim[1]
      maxX <- object@g.args$ylim[2]
    }

    origin <- object@adeg.par$porigin
    lim <- setlimits1D(minX, maxX, origin = origin$origin[1], includeOr = origin$include)

    ## compute grid size
    tmp <- pretty(lim, n = object@adeg.par$pgrid$nint)
    if(!origin$include)
      origin$origin[1] <- tmp[1]
    
    cgrid <- diff(tmp)[1]
    if(is.na(cgrid))
      stop("error while calculating grid")

    ## compute grid location
    v0 <- origin$origin[1]
    if((origin$origin[1] + cgrid) <= lim[2])
      v0 <- c(v0, seq(origin$origin[1] + cgrid, lim[2], by = cgrid))
    if((origin$origin[1] - cgrid >= lim[1]))
      v0 <- c(v0, seq(origin$origin[1] - cgrid, lim[1], by = -cgrid))
    v0 <- sort(v0[v0 >= lim[1] & v0 <= lim[2]])

    ## clean near-zero values
    delta <- diff(range(v0))/object@adeg.par$pgrid$nint
    if (any(small <- abs(v0) < 1e-14 * delta)) 
        v0[small] <- 0
    
    object@s.misc$backgrid <- list(x = v0, d = cgrid)
    
    ## object@adeg.par$paxes has priority over object@g.args$scales
    scalesandlab <- modifyList(as.list(object@g.args$scales), object@adeg.par$paxes, keep.null = TRUE)
    
    if(!scalesandlab$draw) {
      scalesandlab$x$draw <- FALSE
      scalesandlab$y$draw <- FALSE
    }
    
    lead <- ifelse(object@adeg.par$p1d$reverse, 1 , -1)
    
    if(object@adeg.par$p1d$horizontal) {
      ## draw axes for horizontal plot
      if(is.null(scalesandlab$x$at))
        scalesandlab$x$at <- object@s.misc$backgrid$x
      
      if(is.null(object@g.args$xlim) || !identical(object@s.misc$hori.update, object@adeg.par$p1d$horizontal))
        object@g.args$xlim <- lim
      
      if(is.null(object@g.args$ylim))
        object@g.args$ylim <- setlimits1D(min(at), max(at), 0, FALSE)
      if(inherits(object, "S1.boxplot")) ## extend ylim for boxes
        object@g.args$ylim <- object@g.args$ylim + c(-1, 1) * abs(diff(range(at))) / (nlevels(fac) + 1)
      
      ref <- ifelse(object@adeg.par$p1d$reverse, 2, 1)
      margin <- object@g.args$ylim[ref]
      if(object@adeg.par$p1d$rug$draw)
        margin <- object@adeg.par$p1d$rug$margin * abs(diff(object@g.args$ylim))
      object@s.misc$rug <- object@g.args$ylim[ref]
      object@g.args$ylim[ref] <- object@g.args$ylim[ref] + lead * margin
        
    } else {
      ## draw axes for vertical plot
      if(is.null(scalesandlab$y$at))
        scalesandlab$y$at <- object@s.misc$backgrid$x
      
      if(is.null(object@g.args$ylim) || !identical(object@s.misc$hori.update, object@adeg.par$p1d$horizontal))
        object@g.args$ylim <- lim
      
      if(is.null(object@g.args$xlim))
        object@g.args$xlim <- setlimits1D(min(at), max(at), 0, FALSE)
      if(inherits(object, "S1.boxplot")) ## extend xlim for boxes
        object@g.args$xlim <- object@g.args$xlim + c(-1, 1) * abs(diff(range(at))) / (nlevels(fac) + 1)
      
      ref <- ifelse(object@adeg.par$p1d$reverse, 2, 1)
      margin <- object@g.args$xlim[ref]
      if(object@adeg.par$p1d$rug$draw)
        margin <- object@adeg.par$p1d$rug$margin * abs(diff(object@g.args$xlim))
      object@s.misc$rug <- object@g.args$xlim[ref]
      object@g.args$xlim[ref] <-  object@g.args$xlim[ref] + lead * margin
    }
    
    object@g.args$scales <- scalesandlab
    object@s.misc$hori.update <- object@adeg.par$p1d$horizontal
    assign(name_obj, object, envir = parent.frame())
  })


setMethod(
  f = "panelbase",
  signature = "ADEg.S1",
  definition = function(object, x, y) {
    ## Formula defined in gettrellis
    ## if horizontal, x is score and y is a vector with repetitions of origin
    ## if vertical, this is the inverse
    
    grid <- object@adeg.par$pgrid
    porigin <- object@adeg.par$porigin 
    pscore <- object@adeg.par$p1d
    lims <- current.panel.limits(unit = "native")

    plines <- object@adeg.par$plines
    if(!is.null(object@data$fac)) {
      ## there is a factor in the data (e.g., S1.class)
      if(object@data$storeData)
        fac <- object@data$fac
      else
        fac <- as.factor(eval(object@data$fac, envir = sys.frame(object@data$frame)))
      plines <- lapply(plines, FUN = function(x) return(rep(x, length.out = nlevels(fac))[fac]))
    }
    
    lead <- ifelse(pscore$reverse, -1 , 1)
    
    if(pscore$horizontal) {
      ## horizontal plot
      
      ## set margins to get some place for rug
      ref <- ifelse(pscore$reverse, object@g.args$ylim[2], object@g.args$ylim[1])
      margin <- ref
      if(pscore$rug$draw)
        margin <- ifelse(is.unit(pscore$rug$margin), convertUnit(pscore$rug$margin, typeFrom = "dimension", unitTo = "native", axisFrom = "y", valueOnly = TRUE), pscore$rug$margin)
            
      ## draw grid
      if(grid$draw)
        panel.segments(x0 = object@s.misc$backgrid$x , x1 = object@s.misc$backgrid$x, y0 = lims$ylim[1], y1 = lims$ylim[2], col = grid$col, lty = grid$lty, lwd = grid$lwd)
      
      ## draw origin
      panel.abline(
        v = if(porigin$draw) porigin$origin else NULL,
        h = if(pscore$rug$draw & pscore$rug$line) object@s.misc$rug else NULL,
        col = porigin$col, lwd = porigin$lwd, lty = porigin$lty, alpha = porigin$alpha)
      
      ## draw rug
      if(pscore$rug$draw & (pscore$rug$tck != 0)) {
        ## tick end and starting points
        start <- object@s.misc$rug
        end <- start - pscore$rug$tck * lead * abs(start - ref)
        ## 'panel.rug' needs 'npc' values 
        start <- convertUnit(unit(start, "native"), unitTo = "npc", axisFrom = "y", valueOnly = TRUE)
        end <- convertUnit(unit(end, "native"), unitTo = "npc", axisFrom = "y", valueOnly = TRUE)
        do.call("panel.rug", c(list(x = y, start = start, end = end), plines))
      }
    } else {
      ## vertical plot
      
      ## set margins to get some place for rug
      ref <- ifelse(pscore$reverse, object@g.args$xlim[2], object@g.args$xlim[1])
      margin <- ref
      if(pscore$rug$draw)          
        margin <- ifelse(is.unit(pscore$rug$margin), convertUnit(pscore$rug$margin, typeFrom = "dimension", unitTo = "native", axisFrom = "x", valueOnly = TRUE), pscore$rug$margin)
      
      ## draw grid
      if(grid$draw)
        panel.segments(y0 = object@s.misc$backgrid$x , y1 = object@s.misc$backgrid$x, x0 = lims$xlim[1], x1 = lims$xlim[2], col = grid$col, lty = grid$lty, lwd = grid$lwd)

      ## draw origin
      panel.abline(
        h = if(porigin$draw) porigin$origin else NULL,
        v = if(pscore$rug$draw & pscore$rug$line) object@s.misc$rug else NULL,
        col = porigin$col, lwd = porigin$lwd, lty = porigin$lty, alpha = porigin$alpha)

      ## draw rug
      if(pscore$rug$draw && (pscore$rug$tck != 0)) {
        ## tick end and starting points
        start <- object@s.misc$rug
        end <- start - pscore$rug$tck * lead * abs(start - ref)
        start <- convertUnit(unit(start, "native"), unitTo = "npc", axisFrom = "x", valueOnly = TRUE)
        end <- convertUnit(unit(end, "native"), unitTo = "npc", axisFrom = "x", valueOnly = TRUE)
        do.call("panel.rug", c(list(y = y, start = start, end = end), plines))
      }
    }

    ## indicate grid size (d = **)
    if(grid$draw & (grid$text$cex > 0)) { 
      text.pos <- .setposition(grid$text$pos)
      textgrid <- textGrob(label = paste("d =", object@s.misc$backgrid$d), x = text.pos$posi[1], y = text.pos$posi[2], gp = gpar(cex = grid$text$cex, col = grid$text$col), name = "gridtext")
      grid.rect(x = text.pos$posi[1], y = text.pos$posi[2], width = grobWidth(textgrid), height = grobHeight(textgrid), gp = gpar(fill= object@adeg.par$pbackground$col, alpha = 0.8, col = "transparent"))
      grid.draw(textgrid)
    }
    callNextMethod()
  })


setMethod(
  f = "setlatticecall",
  signature = "ADEg.S1",
  definition = function(object) {
    ## arguments recurrents de la liste, pas les limites car elles seront definis ensuite
    name_obj <- deparse(substitute(object))
    
    ## grid background and box
    object@trellis.par$panel.background$col <- object@adeg.par$pbackground$col
    if(!object@adeg.par$pbackground$box)
      object@trellis.par$axis.line$col <- "transparent"
    
    arguments <- list(
                   par.settings = object@trellis.par,
                   scales = object@g.args$scales,
                   key = createkey(object),                   
                   axis = axis.L, ## see utils.R
                   panel = function(...) {
                     panelbase(object,...) ## grid,
                     panel(object,...) ## call to S1.panel function, for slabel and ADEg.S1 class of graphs
                   })
    
    object@lattice.call$arguments <- arguments
    object@lattice.call$graphictype <- "xyplot" 

    ## get lattice arguments (set unspecified to NULL)
    argnames <- c("main", "sub", "xlab", "ylab")
    largs <- object@g.args[argnames]
    names(largs) <- argnames
    ## add xlim and ylim if not NULL
    if("xlim" %in% names(object@g.args))
      largs["xlim"] <- object@g.args["xlim"]
    if("ylim" %in% names(object@g.args))
      largs["ylim"] <- object@g.args["ylim"]
    
    object@lattice.call$arguments <- c(object@lattice.call$arguments, largs, list(strip = FALSE))
    assign(name_obj, object, envir = parent.frame())
})


setMethod(
  f = "gettrellis",
  signature = "ADEg.S1",
  definition = function(object) {
    if(object@data$storeData)
      score <- object@data$score
    else
      score <- eval(object@data$score, envir = sys.frame(object@data$frame))
    
    score <- as.matrix(score)[, 1]  ## to manage 'score' when it is a data.frame with only one column
    
    xdata <- rep(1, length(score))
    fml <- as.formula(score ~ xdata)
    tmptrellis <- do.call(what = object@lattice.call$graphictype, args = c(fml, object@lattice.call$arguments, environment()))
    return(tmptrellis)
  })


## zoom without center
setMethod(
  f = "zoom",
  signature = c("ADEg.S1", "numeric", "missing"),
  definition = function(object, zoom, center) {
    ## zoom in xlim
    p1d <- object@adeg.par$p1d
    nameobj <- deparse(substitute(object))
    if(length(zoom) != 1)
      stop("zoom factor should be length 1")
    center <- ifelse(p1d$horizontal, mean(object@g.args$xlim), mean(object@g.args$ylim))
    zoom(object, zoom, center)
  })


## zoom with center
setMethod(
  f = "zoom",
  signature = c("ADEg.S1",  "numeric", "numeric"),
  definition = function(object, zoom, center) {
    nameobj <- deparse(substitute(object))
    p1d <- object@adeg.par$p1d
    origin <- object@adeg.par$porigin
    
    if(length(center) != 1)
      stop("Center should be a numeric")
    if(length(zoom) != 1)
      stop("Zoom factor should be a numeric")
    
    if(p1d$horizontal) {
      diffx <- diff(object@g.args$xlim) / zoom
      minX <- center - diffx / 2
      maxX <- center + diffx / 2
      object@g.args$xlim <- c(minX, maxX)
    } else {
      diffx <- diff(object@g.args$ylim) / zoom
      minX <- center - diffx / 2
      maxX <- center + diffx / 2
      object@g.args$ylim <- c(minX, maxX)
    }
    
    lim <- setlimits1D(minX, maxX, origin = origin$origin[1], includeOr = origin$include)
    ## compute grid size
    tmp <- pretty(lim, n = object@adeg.par$pgrid$nint)
    if(!origin$include)
      origin$origin[1] <- tmp[1]
    
    cgrid <- diff(tmp)[1]
    if(is.na(cgrid))
      stop("error while calculating grid")
    
    ## compute grid location
    v0 <- origin$origin[1]
    if((origin$origin[1] + cgrid) <= lim[2])
      v0 <- c(v0, seq(origin$origin[1] + cgrid, lim[2], by = cgrid))
    if((origin$origin[1] - cgrid >= lim[1]))
      v0 <- c(v0, seq(origin$origin[1] - cgrid, lim[1], by = -cgrid))
    v0 <- sort(v0[v0 >= lim[1] & v0 <= lim[2]])
    object@s.misc$backgrid <- list(x = v0, d = cgrid)

    setlatticecall(object)
    print(object)
    invisible()   
  })
