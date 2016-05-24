####################################################
##              Curves Plot                       ##
##        1d score represents in 2D plot          ##
####################################################
setClass(
    Class = "ADEg.C1",
    contains = c("ADEg", "VIRTUAL"),
    slots = c(data = "list")
    )

setMethod(
    f = "initialize",
    signature  = "ADEg.C1",
    definition = function(.Object, data = list(score = NULL, frame = 0, storeData = TRUE), ...) {
        .Object <- callNextMethod(.Object, ...) ## ADEg initialize
        .Object@data <- data
        return(.Object)
    })


setMethod(
    f = "prepare",
    signature = "ADEg.C1",
    definition = function(object) {
        ## prepare: grid calculations
        ## reset limits and sets axis information for lattice
        
        name_obj <- deparse(substitute(object))
        if(object@data$storeData)
            score <- object@data$score
        else
            score <- eval(object@data$score, envir = sys.frame(object@data$frame))
        
        if(inherits(object, "C1.curve") | inherits(object, "C1.dotplot") | inherits(object, "C1.interval"))
            if(object@data$storeData)
                at <- object@data$at
            else
                at <- eval(object@data$at, envir = sys.frame(object@data$frame))
        
        score <- as.matrix(score)[, 1]  ## to manage 'score' when it is a data.frame with only one column
        
        if(class(object) == "C1.interval")  ## to manage only the first score in c(score1, score2)
            score <- score[1:(length(score) / 2)]
        
        ## limits and scale
        if(!is.null(object@s.misc$hori.update))
            if(object@s.misc$hori.update != object@adeg.par$p1d$horizontal) {
                aux <- object@g.args$xlim
                object@g.args$xlim <- object@g.args$ylim
                object@g.args$ylim <- aux 
            }
        object@s.misc$hori.update <- object@adeg.par$p1d$horizontal
        
        minX <- min(score)
        maxX <- max(score)
        if(object@adeg.par$p1d$horizontal & !is.null(object@g.args$xlim)) {
            minX <- object@g.args$xlim[1]
            maxX <- object@g.args$xlim[2]
        }

        if(!object@adeg.par$p1d$horizontal & !is.null(object@g.args$ylim)) {
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
            
            if(is.null(object@g.args$xlim))
                object@g.args$xlim <- lim
            
            if(inherits(object, "C1.curve") | inherits(object, "C1.dotplot") | inherits(object, "C1.interval"))
                if(!is.null(at))
                    scalesandlab$y$at <- at
            
            if(is.null(scalesandlab$y$at))
                scalesandlab$y$at <- 1:NROW(score)
            
        } else {
            ## draw axes for vertical plot
            if(is.null(scalesandlab$y$at))
                scalesandlab$y$at <- object@s.misc$backgrid$x
            
            if(is.null(object@g.args$ylim))
                object@g.args$ylim <- lim
            
            if(inherits(object, "C1.curve") | inherits(object, "C1.dotplot") | inherits(object, "C1.interval"))
                if(!is.null(at))
                    scalesandlab$x$at <- at
            
            if(is.null(scalesandlab$x$at))
                scalesandlab$x$at <- 1:NROW(score)
        }
        
        object@g.args$scales <- scalesandlab
        assign(name_obj, object, envir = parent.frame())
    })


setMethod(
    f = "panelbase",
    signature = "ADEg.C1",
    definition = function(object, x, y) {
        ## Formula defined in gettrellis
        ## if horizontal, x is score and y is a vector with repetitions of origin
        ## if vertical, this is the inverse
        grid <- object@adeg.par$pgrid
        porigin <- object@adeg.par$porigin 
        pscore <- object@adeg.par$p1d
        lims <- current.panel.limits(unit = "native")
        
        ## for rugs
        if(pscore$rug$draw & (pscore$rug$tck != 0)) {
            plines <- object@adeg.par$plines
            if(!is.null(object@data$fac)) {
                ## C1.density or C1.gauss (different colors for rugs)
                if(object@data$storeData)
                    fac <- as.factor(object@data$fac)
                else
                    fac <- as.factor(eval(object@data$fac, envir = sys.frame(object@data$frame)))
                plines <- lapply(plines, FUN = function(x) {return(rep(x, length.out = nlevels(fac))[fac])})
            }
        }
        lead <- ifelse(pscore$reverse, -1, 1)
        
        if(pscore$horizontal) {
            ## horizontal plot
            
            ## set margins to get some place for rug
            ref <- ifelse(pscore$reverse, lims$ylim[2], lims$ylim[1])
            margin <- ref
            if(pscore$rug$draw)
                margin <- ifelse(is.unit(pscore$rug$margin), convertUnit(pscore$rug$margin, typeFrom = "dimension", unitTo = "native", axisFrom = "y", valueOnly = TRUE), pscore$rug$margin)
            
            ## draw grid
            if(grid$draw)
                panel.segments(x0 = object@s.misc$backgrid$x , x1 = object@s.misc$backgrid$x, y0 = lims$ylim[1], y1 = lims$ylim[2], col = grid$col, lty = grid$lty, lwd = grid$lwd)
            
            ## draw origin
            panel.abline(
                v = if(porigin$draw) porigin$origin else NULL,
                h = if(pscore$rug$draw & pscore$rug$line) ref + lead * margin else NULL,
                col = porigin$col, lwd = porigin$lwd, lty = porigin$lty, alpha = porigin$alpha)
            
            ## draw rug
            if(pscore$rug$draw & (pscore$rug$tck != 0)) {
                ## tick end and starting points
                start <- ref + lead * margin
                end <- (start - pscore$rug$tck * lead * abs(start - ref))
                start <- convertUnit(unit(start, "native"), unitTo = "npc", axisFrom = "y", valueOnly = TRUE)
                end <- convertUnit(unit(end, "native"), unitTo = "npc", axisFrom = "y", valueOnly = TRUE)
                do.call("panel.rug", c(list(x = y, start = start, end = end), plines))
            }
            
        } else {
            ## vertical plot
            
            ## set margins to get some place for rug
            ref <- ifelse(pscore$reverse, lims$xlim[2], lims$xlim[1])
            margin <- ref
            if(pscore$rug$draw)          
                margin <- ifelse(is.unit(pscore$rug$margin), convertUnit(pscore$rug$margin, typeFrom = "dimension", unitTo = "native", axisFrom = "x", valueOnly = TRUE), pscore$rug$margin)
            
            ## draw grid
            if(grid$draw)
                panel.segments(y0 = object@s.misc$backgrid$x , y1 = object@s.misc$backgrid$x, x0 = lims$xlim[1], x1 = lims$xlim[2], col = grid$col, lty = grid$lty, lwd = grid$lwd)

            ## draw origin
            panel.abline(
                h = if(porigin$draw) porigin$origin else NULL,
                v = if(pscore$rug$draw & pscore$rug$line) ref + lead * margin else NULL,
                col = porigin$col, lwd = porigin$lwd, lty = porigin$lty, alpha = porigin$alpha)

            ## draw rug
            if(pscore$rug$draw && pscore$rug$tck != 0) {
                ## tick end and starting points
                start <- ref + lead * margin
                end <- (start - pscore$rug$tck * lead * abs(start - ref))
                start <- convertUnit(unit(start, "native"), unitTo = "npc", axisFrom = "x", valueOnly = TRUE)
                end <- convertUnit(unit(end, "native"), unitTo = "npc", axisFrom = "x", valueOnly = TRUE)
                do.call("panel.rug", c(list(y = y, start = start, end = end), plines))
            }
        }

        ## indicate grid size (d = **)
        if(grid$draw & (grid$text$cex > 0)) { 
            text.pos <- .setposition(grid$text$pos)
            textgrid <- textGrob(label = paste("d =", object@s.misc$backgrid$d), x = text.pos$posi[1], y = text.pos$posi[2], gp = gpar(cex = grid$text$cex, col = grid$text$col), name = "gridtext")
            grid.rect(x = text.pos$posi[1], y = text.pos$posi[2], width = grobWidth(textgrid), height = grobHeight(textgrid), gp = gpar(fill = object@adeg.par$pbackground$col, alpha = 0.8, col = "transparent"))
            grid.draw(textgrid)
        }
        
        callNextMethod()
    })


setMethod(
    f = "setlatticecall",
    signature = "ADEg.C1",
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
                panelbase(object, ...) ## grid,
                panel(object, ...) ## call to C1.panel function, for slabel and ADEg.C1 class of graphs
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
    signature = "ADEg.C1",
    definition = function(object) {
        if(object@data$storeData)
            score <- object@data$score
        else
            score <- eval(object@data$score, envir = sys.frame(object@data$frame))
        
        score <- as.matrix(score)[, 1]  ## to manage 'score' when it is a data.frame with only one column

        xdata <- rep(1, length(score))
        if(inherits(object, "C1.barchart")) {
            xdata <- 1:length(score)
        } else if(inherits(object, "C1.dotplot") | inherits(object, "C1.curve") | inherits(object, "C1.interval")) {
            if(object@data$storeData)
                xdata <- object@data$at
            else
                xdata <- eval(object@data$at, envir = sys.frame(object@data$frame))
        }
        
        fml <- as.formula(score ~ xdata)
        
        tmptrellis <- do.call(what = object@lattice.call$graphictype, args = c(fml, object@lattice.call$arguments, environment()))
        return(tmptrellis)
    })
