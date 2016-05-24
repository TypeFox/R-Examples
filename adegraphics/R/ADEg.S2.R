setClass(
    Class = "ADEg.S2",
    contains = c("ADEg", "VIRTUAL"),
    slots = c(data = "list")
    )


setMethod(
    f = "initialize",
    signature = "ADEg.S2",
    definition = function(.Object, data = list(dfxy = NULL, xax = 1, yax = 2, frame = 0, storeData = TRUE), ...) {
        .Object <- callNextMethod(.Object, ...) ## ADEg initialize
        .Object@data <- data
        return(.Object)
    })


setMethod(
    f = "prepare",
    signature = "ADEg.S2",
    definition = function(object) {
        ## TODO: factorise les if
        name_obj <- deparse(substitute(object))
        if(object@data$storeData)
            dfxy <- object@data$dfxy
        else
            dfxy <- eval(object@data$dfxy, envir = sys.frame(object@data$frame))
        
        ## axes limits
        if(class(object) == "S2.corcircle") {
            object@trellis.par$panel.background$col <- "transparent"
            if(object@g.args$fullcircle) {
                if(is.null(object@g.args$xlim) || !identical(object@s.misc$xfullcircle.update, object@g.args$fullcircle)) {
                    minX <- -1
                    maxX <- 1
                } else {
                    minX <- object@g.args$xlim[1]
                    maxX <- object@g.args$xlim[2]
                }
                if(is.null(object@g.args$ylim) || !identical(object@s.misc$yfullcircle.update, object@g.args$fullcircle)) {
                    minY <- -1
                    maxY <- 1
                } else {
                    minY <- object@g.args$ylim[1]
                    maxY <- object@g.args$ylim[2]
                }
            } else {
                if(is.null(object@g.args$xlim) || !identical(object@s.misc$xfullcircle.update, object@g.args$fullcircle)) {
                    minX <- min(dfxy[, object@data$xax])
                    maxX <- max(dfxy[, object@data$xax])
                } else {
                    minX <- object@g.args$xlim[1]
                    maxX <- object@g.args$xlim[2]
                }
                if(is.null(object@g.args$ylim) || !identical(object@s.misc$yfullcircle.update, object@g.args$fullcircle)) {
                    minY <- min(dfxy[, object@data$yax])
                    maxY <- max(dfxy[, object@data$yax])
                } else {
                    minY <- object@g.args$ylim[1]
                    maxY <- object@g.args$ylim[2]
                }
            }
        } else {
            if(is.null(object@g.args$xlim)) {
                minX <- min(dfxy[, object@data$xax])
                maxX <- max(dfxy[, object@data$xax])
            } else {
                minX <- object@g.args$xlim[1]
                maxX <- object@g.args$xlim[2]
            }
            if(is.null(object@g.args$ylim)) {
                minY <- min(dfxy[, object@data$yax])
                maxY <- max(dfxy[, object@data$yax])
            } else {
                minY <- object@g.args$ylim[1]
                maxY <- object@g.args$ylim[2]
            }
        }
        
        limits <- setlimits2D(minX = minX, maxX = maxX, minY = minY, maxY = maxY, origin = rep(object@adeg.par$porigin$origin, le = 2),
                             aspect.ratio = object@adeg.par$paxes$aspectratio, includeOr = object@adeg.par$porigin$include)
        
        if(is.null(object@g.args$xlim) || !identical(object@s.misc$xfullcircle.update, object@g.args$fullcircle))
            object@g.args$xlim <- limits$xlim
        if(is.null(object@g.args$ylim) || !identical(object@s.misc$yfullcircle.update, object@g.args$fullcircle))
            object@g.args$ylim <- limits$ylim
        
        if(class(object) == "S2.corcircle") {
            object@s.misc$xfullcircle.update <- object@g.args$fullcircle
            object@s.misc$yfullcircle.update <- object@g.args$fullcircle
        }
        
        ## grid locations and axes 
        if(object@adeg.par$pgrid$draw || object@adeg.par$paxes$draw) {
            ## axes division
            if(class(object) != "S2.corcircle") {
                if(object@adeg.par$porigin$include)
                    object@s.misc$backgrid <- .getgrid(xlim = object@g.args$xlim, ylim = object@g.args$ylim, object@adeg.par$pgrid$nint, rep(object@adeg.par$porigin$origin, le = 2), asp = object@adeg.par$paxes$aspectratio)
                else
                    object@s.misc$backgrid <- .getgrid(xlim = object@g.args$xlim, ylim = object@g.args$ylim, object@adeg.par$pgrid$nint, asp = object@adeg.par$paxes$aspectratio)
            }
            
            if(object@adeg.par$paxes$draw) {
                ## parameters to plot axes
                scalesandlab <- list(x = object@adeg.par$paxes$x, y = object@adeg.par$paxes$y)
                if(is.null(scalesandlab$y$at)) {
                    scalesandlab$y$at <- object@s.misc$backgrid[[3L]][!is.na(object@s.misc$backgrid[[3L]])]
                    if(class(object) == "S2.corcircle")
                        scalesandlab$y$at <- scalesandlab$y$at[(length(scalesandlab$y$at) / 2 + 1):length(scalesandlab$y$at)]
                }
                if(is.null(scalesandlab$x$at)) {
                    scalesandlab$x$at <- object@s.misc$backgrid[[1L]][!is.na(object@s.misc$backgrid[[1L]])]
                    if(class(object) == "S2.corcircle")
                        scalesandlab$x$at <- scalesandlab$x$at[1:(length(scalesandlab$x$at) / 2)]
                }
            } else 
                scalesandlab <- list(draw = FALSE) ## no axes
        }
        else
            scalesandlab <- list(draw = FALSE) ## no axes
        
        if(object@adeg.par$paxes$aspectratio != "iso")
            object@adeg.par$pgrid$text$cex <- 0 ## grid cell size has no meaning
        
        ## if grid and axes are drawn, no text indication
        if(object@adeg.par$pgrid$draw && object@adeg.par$paxes$draw)
            object@adeg.par$pgrid$text$cex <- 0
        object@g.args$scales <- scalesandlab
        assign(name_obj, object, envir = parent.frame())
    })


setMethod(
    f = "panelbase",
    signature = "ADEg.S2",
    definition = function(object, x, y) {
        ## draw grid
        lims <- current.panel.limits(unit = "native")
        porigin <- object@adeg.par$porigin
        porigin$origin <- rep(porigin$origin, length.out = 2)

        if(class(object) == "S2.corcircle") 
            grid.circle(x = 0, y = 0, r = 1, default.units = "native", gp = gpar(col = "black", fill = object@adeg.par$pbackground$col), draw = TRUE, name = "circleGrid")
        
        if(object@adeg.par$pgrid$draw) { ## if grid to draw
            grid <- object@adeg.par$pgrid
            locations <- object@s.misc$backgrid ## coordinates for the grid 
            panel.segments(
                x0 = c(locations$x0[!is.na(locations$x0)], rep(lims$xlim[1], sum(is.na(locations$x0)))),
                x1 = c(locations$x1[!is.na(locations$x1)], rep(lims$xlim[2], sum(is.na(locations$x1)))),
                y0 = c(rep(lims$ylim[1], sum(is.na(locations$y0))), locations$y0[!is.na(locations$y0)]),
                y1 = c(rep(lims$ylim[2], sum(is.na(locations$y1))), locations$y1[!is.na(locations$y1)]),
                col = grid$col, lty = grid$lty, lwd = grid$lwd)
            
            if(grid$text$cex > 0) {
                text.pos <- .setposition(grid$text$pos)
                textgrid <- textGrob(label = paste("d =", locations$d), x = text.pos$posi[1], y = text.pos$posi[2], just = text.pos$just, gp = gpar(cex = grid$text$cex, col = grid$text$col), name = "gridtext")
                grid.rect(x = text.pos$posi[1], y = text.pos$posi[2], width = grobWidth(textgrid), height = grobHeight(textgrid),
                          just = text.pos$just, gp = gpar(fill = ifelse(class(object) == "S2.corcircle", "transparent", object@adeg.par$pbackground$col), alpha = 1, col = "transparent"))
                grid.draw(textgrid)
            }
        }
        
        if(porigin$draw && porigin$include & class(object) == "S2.corcircle") {
            panel.segments(x0 = c(-1, porigin$origin[1]), x1 = c(1, porigin$origin[1]), y0 = c(porigin$origin[2], -1), y1 = c(porigin$origin[2], 1), col = porigin$col, lwd = porigin$lwd, lty = porigin$lty, alpha = porigin$alpha)
            ## TODO: check last parameters valididy     
        }
        
        if(porigin$draw && porigin$include & !class(object) == "S2.corcircle") {
            panel.abline(h = porigin$origin[2], v = porigin$origin[1], col = porigin$col, lwd = porigin$lwd, lty = porigin$lty, alpha = porigin$alpha)
            ## TODO: check last parameters valididy
        }
        
        ## spatial object management
        if(any(names(object@g.args) == "Sp")) {
            do.call("adeg.panel.Spatial", args = c(list(SpObject = object@g.args$Sp, sp.layout = object@g.args$sp.layout), object@adeg.par$pSp))
        }
        else  ## no Sp but sp.layout
            if(any(names(object@g.args) == "sp.layout"))
              sppanel(lst = object@g.args$sp.layout)
        
        ## neighbouring object management
        if(any(names(object@g.args) == "nbobject")) {
            nbobj <- object@g.args$nbobject
            if(class(nbobj) != "nb")
                stop("wrong class for the nb object")
            pnb <- object@adeg.par$pnb
            do.call("adeg.panel.nb", args = list(nbobject = nbobj, coords = cbind(x, y), col.edge = pnb$edge$col, lwd = pnb$edge$lwd, lty = pnb$edge$lty, pch = pnb$node$pch, cex = pnb$node$cex, col.node = pnb$node$col, alpha = pnb$node$alpha))
        }
        callNextMethod()
    })


setMethod(
    f = "setlatticecall",
    signature = "ADEg.S2",
    definition =  function(object) {
        ## arguments recurrents de la liste, pas les limites car elles seront definis ensuite
        name_obj <- deparse(substitute(object))

        ## background and box
        if(!inherits(object, "S2.corcircle"))
            object@trellis.par$panel.background$col <- object@adeg.par$pbackground$col
        if(!object@adeg.par$pbackground$box)
            object@trellis.par$axis.line$col <- "transparent"
        
        arguments <- list(
            par.settings = object@trellis.par,
            scales = object@g.args$scales,
            aspect = object@adeg.par$paxes$aspectratio,
            key = createkey(object),
            legend = createcolorkey(object),
            axis = axis.L, ## see utils.R
            panel = function(...) {
                panelbase(object,...) ## grid,
                panel(object,...) ## call to S2.panel function, for slabel and ADEg.S2 class of graphs
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


## zoom without center
setMethod(
    f = "zoom",
    signature = c("ADEg.S2", "numeric", "missing"),
    definition  =  function(object, zoom, center) {
        oldxlim <- object@g.args$xlim
        oldylim <- object@g.args$ylim
        if(length(zoom) != 1)
            stop("zoom factor should be length 1")
        diffx <- diff(oldxlim)
        diffy <- diff(oldylim)
        center <- c(oldxlim[1] + diffx / 2, oldylim[1] + diffy / 2)
        diffx <- diffx / zoom
        diffy <- diffy / zoom
        object@g.args$xlim <- c(center[1] - diffx / 2, center[1] + diffx / 2)
        object@g.args$ylim <- c(center[2] - diffy / 2, center[2] + diffy / 2)
        if(object@adeg.par$pgrid$draw || object@adeg.par$paxes$draw)
            object@s.misc$backgrid <- .getgrid(xlim = object@g.args$xlim, ylim = object@g.args$ylim, object@adeg.par$pgrid$nint, object@adeg.par$porigin$origin, asp = object@adeg.par$paxes$aspectratio)
        prepare(object)
        setlatticecall(object)
        print(object)
        invisible()
    })


## zoom with center
setMethod(
    f = "zoom",
    signature = c("ADEg.S2", "numeric", "numeric"),
    definition = function(object, zoom, center) {
        if(length(center) != 2) 
            stop("error, center should be length 2")
        if(length(zoom) != 1) 
            stop("zoom factor should be length 1")
        diffx <- diff(object@g.args$xlim) / zoom
        diffy <- diff(object@g.args$ylim) / zoom
        object@g.args$xlim <- c(center[1] - diffx / 2, center[1] + diffx / 2)
        object@g.args$ylim <- c(center[2] - diffy / 2, center[2] + diffy / 2)
        if(object@adeg.par$pgrid$draw || object@adeg.par$paxes$draw)
            object@s.misc$backgrid <- .getgrid(xlim = object@g.args$xlim, ylim = object@g.args$ylim, object@adeg.par$pgrid$nint, object@adeg.par$porigin$origin, asp = object@adeg.par$paxes$aspectratio)
        prepare(object)
        setlatticecall(object)
        print(object)
        invisible()
    })


setMethod(
    f = "addhist",
    signature = "ADEg.S2",
    definition = function(object, bandwidth, gridsize = 60, kernel = "normal", cbreaks = 2, storeData = TRUE, plot = TRUE, pos = -1, ...) {
        thecall <- .expand.call(match.call())
        dfcall <- thecall$object
        dfxycall <- substitute(dfcall@data$dfxy)
        
        if(!(inherits(object, "ADEg.S2")))
            stop("Only implemented for 'ADEg.S2' object")
        
        if(storeData) {
            dfxy <- object@data$dfxy
            xax <- object@data$xax
            yax <- object@data$yax
        } else {
            dfxy <- eval(object@data$dfxy, envir = sys.frame(object@data$frame))
            xax <- eval(object@data$xax, envir = sys.frame(object@data$frame))
            yax <- eval(object@data$yax, envir = sys.frame(object@data$frame))
        }
        
        ## sorting parameters
        graphsnames <- c("object", "densX", "densY", "link") 
        sortparameters <- sortparamADEgS(..., graphsnames = graphsnames)
        params <- vector("list", 4)
        names(params) <- graphsnames
        sortparameters <- modifyList(params, sortparameters, keep.null = TRUE)
        update(object, sortparameters[[1]], plot = FALSE)
        
        ## setting positions
        positions <- layout2position(matrix(c(2, 4, 1, 3), 2, 2, byrow = TRUE), c(3, 1) / 2, c(3, 1) / 2, FALSE)

        ## grid computation
        xlimX <- object@g.args$xlim
        ylimY <- object@g.args$ylim
        breaks <- object@s.misc$backgrid
        cgrid <- breaks$d / cbreaks
        bb1 <- range(breaks$x0[!is.na(breaks$x0)])
        bb2 <- range(breaks$y0[!is.na(breaks$y0)])
        breaksX <- seq(from = bb1[1], to = bb1[2], by = cgrid)
        breaksY <- seq(from = bb2[1], to = bb2[2], by = cgrid)
        while(min(breaksX) > xlimX[1])
            breaksX <- c((min(breaksX) - cgrid), breaksX)
        while(max(breaksX) < xlimX[2])
            breaksX <- c(breaksX, (max(breaksX) + cgrid))
        while(min(breaksY) > ylimY[1])
            breaksY <- c((min(breaksY) - cgrid), breaksY)
        while(max(breaksY) < ylimY[2])
            breaksY <- c(breaksY, (max(breaksY) + cgrid))
        
        ## limits and graduation
        dfxaxcall <- call("[", dfxycall, 1:NROW(eval(dfxycall)), substitute(xax))
        dfxaxcallplus <- call("~", dfxaxcall, 1)
        dfyaxcall <- call("[", dfxycall, 1:NROW(eval(dfxycall)), substitute(yax))
        dfyaxcallplus <- call("~", dfyaxcall, 1)
        limcalX <- hist(dfxy[, xax], breaksX, plot = FALSE)
        limcalXcall <- call("hist", substitute(dfxaxcall), breaksX, plot = FALSE)
        limcalY <- hist(dfxy[, yax], breaksY, plot = FALSE)
        limcalYcall <- call("hist", substitute(dfyaxcall), breaksY, plot = FALSE)
        
        top <- 1.1 * max(c(limcalX$counts, limcalY$counts))
        xlimY <- ylimX <- c(0, top)
        drawLines <- pretty(0:top)
        drawLines <- drawLines[-c(1, length(drawLines))]

        if(!missing(bandwidth)) {
            densiX <- bkde(dfxy[, xax], kernel = kernel, bandwidth = bandwidth, gridsize = gridsize)
            densiXcall <- call("bkde", substitute(dfxaxcall), kernel = kernel, bandwidth = bandwidth, gridsize = gridsize)
            densiY <- bkde(dfxy[, yax], kernel = kernel, bandwidth = bandwidth, gridsize = gridsize)
            densiYcall <- call("bkde", substitute(dfyaxcall), kernel = kernel, bandwidth = bandwidth, gridsize = gridsize)
        } else {
            densiX <- bkde(dfxy[, xax], kernel = kernel, gridsize = gridsize)
            densiXcall <- call("bkde", substitute(dfxaxcall), kernel = kernel, gridsize = gridsize)
            densiY <- bkde(dfxy[, yax], kernel = kernel, gridsize = gridsize)
            densiYcall <- call("bkde", substitute(dfyaxcall), kernel = kernel, gridsize = gridsize)
        }
        
        ## trellis creation 
        g2 <- xyplot(dfxy[, xax] ~ 1, xlim = xlimX, ylim = ylimX, horizontal = TRUE, scales = list(draw = FALSE), xlab = NULL, ylab = NULL, histValues = limcalX, 
                     drawLines = drawLines, densi = densiX, params = sortparameters[[2]], 
                     panel = function(histValues, horizontal, drawLines, densi, params) adeg.panel.hist(histValues = histValues, horizontal = horizontal, 
                                                                                         drawLines = drawLines, densi = densi, params = params))
        g2$call <- call("xyplot", dfxaxcallplus, xlim = substitute(xlimX), ylim = substitute(ylimX), horizontal = TRUE, scales = list(draw = FALSE), xlab = NULL, ylab = NULL, 
                        histValues = limcalXcall, drawLines = substitute(drawLines), densi = substitute(densiXcall), params = sortparameters[[2]], 
                        panel = function(histValues, horizontal, drawLines, densi, params) adeg.panel.hist(histValues = histValues, horizontal = horizontal, 
                                                                                            drawLines = drawLines, densi = densi, params = params))
        
        
        g3 <- xyplot(dfxy[, yax] ~ 1, xlim = xlimY, ylim = ylimY, horizontal = FALSE, scales = list(draw = FALSE), xlab = NULL, ylab = NULL, histValues = limcalY, 
                     drawLines = drawLines, densi = densiY, params = sortparameters[[3]], 
                     panel = function(histValues, horizontal, drawLines, densi, params) adeg.panel.hist(histValues = histValues, horizontal = horizontal, 
                                                                                         drawLines = drawLines, densi = densi, params = params))
        g3$call <- call("xyplot", dfyaxcallplus, xlim = substitute(xlimY), ylim = substitute(ylimY), horizontal = FALSE, scales = list(draw = FALSE), xlab = NULL, ylab = NULL, 
                        histValues = limcalYcall, drawLines = substitute(drawLines), densi = substitute(densiYcall), params = sortparameters[[3]], 
                        panel = function(histValues, horizontal, drawLines, densi, params) adeg.panel.hist(histValues = histValues, horizontal = horizontal, 
                                                                                            drawLines = drawLines, densi = densi, params = params))
        
        
        g4 <- xyplot(1 ~ 1, xlim = xlimY, ylim = ylimX, scales = list(draw = FALSE), xlab = NULL, ylab = NULL, drawLines = drawLines, params = sortparameters[[4]], 
                     panel = function(drawLines, params) adeg.panel.join(drawLines = drawLines, params = params))
        g4$call <- call("xyplot", 1 ~ 1, xlim = substitute(xlimY), ylim = substitute(ylimX), scales = list(draw = FALSE), xlab = NULL, ylab = NULL, drawLines = substitute(drawLines), 
                        params = sortparameters[[4]], panel = function(drawLines, params) adeg.panel.join(drawLines = drawLines, params = params))
        
        ## ADEgS creation and display
        obj <- new(Class = "ADEgS", ADEglist = list(object, g2, g3, g4), positions = positions, add = matrix(0, ncol = 4, nrow = 4), Call = match.call())
        names(obj) <- graphsnames
        if(plot)
            print(obj)
        invisible(obj)
    })


setMethod(
    f = "gettrellis",
    signature = "ADEg.S2",
    definition = function(object) {
        if(object@data$storeData) {
            dfxy <- as.matrix(object@data$dfxy)
            xax <- object@data$xax
            yax <- object@data$yax
        } else {
            dfxy <- as.matrix(eval(object@data$dfxy, envir = sys.frame(object@data$frame)))
            yax <- eval(object@data$yax, envir = sys.frame(object@data$frame))
            xax <- eval(object@data$xax, envir = sys.frame(object@data$frame))
        }
        
        tmptrellis <- do.call(what = object@lattice.call$graphictype, args = c(formula(dfxy[, yax] ~ dfxy[, xax]), object@lattice.call$arguments, environment()))
        return(tmptrellis)
    })

