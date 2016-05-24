setClass(
    Class = "ADEg.Tr",
    contains = c("ADEg", "VIRTUAL"),
    slots = c(data = "list")
    )


setMethod(
    f = "initialize",
    signature  = "ADEg.Tr",
    definition = function(.Object, data = list(dfxyz = NULL, frame = 0, storeData = TRUE), ...) {
        .Object <- callNextMethod(.Object, ...) ## ADEg initialize
        .Object@data <- data
        return(.Object)
    })


setMethod(
    f = "prepare",
    signature = "ADEg.Tr",
    definition = function(object) {
        name_obj <- deparse(substitute(object))

        if(object@data$storeData)
            df <- object@data$dfxyz
        else 
            df <- eval(object@data$dfxyz, envir = sys.frame(object@data$frame))

        ## define limits
        if(is.null(object@g.args$xlim))
            object@g.args$xlim <- c(-0.8, 0.8)
        if(is.null(object@g.args$ylim))
            object@g.args$ylim <- c(-0.6, 1)
        
        ## grid computation
        if(is.null(object@g.args$max3d))
            object@g.args$max3d <- .trranges(df = df, adjust = object@g.args$adjust)$maxi
        if(is.null(object@g.args$min3d))
            object@g.args$min3d <- .trranges(df = df, adjust = object@g.args$adjust)$mini

        valuLim <- .trranges(df = df, adjust = object@g.args$adjust, min3 = object@g.args$min3d, max3 = object@g.args$max3d)

        ## coordinates for the triangle vertices
        A <- c(-1 / sqrt(2), -1 / sqrt(6))
        B <- c(1 / sqrt(2), -1 / sqrt(6))
        C <- c(0, 2 / sqrt(6))
        object@s.misc$cornerp <- list(A = A, B = B, C = C)
        
        ## coordinates for grid and axes
        ng <- object@adeg.par$pgrid$nint + 1 ## number of grid lines
        pts1 <- pts2 <- pts3 <- c()
        vdivision <- mapply(FUN = function(min, max) seq(min, max, length.out = ng), min = valuLim$mini, max = valuLim$maxi) ## 3 columns: one per axes    

        ## where to draw the division
        indented <- seq(0, 1, length.out = nrow(vdivision))[-c(1, nrow(vdivision))]

        ## axis 1 (A to B)
        pts1 <- matrix(rep(A, length(indented)), ncol = 2, byrow = TRUE) + indented * (matrix(rep(B, length(indented)), ncol = 2, byrow = TRUE) - matrix(rep(A, length(indented)), ncol = 2, byrow = TRUE)) 
        ##axis 2 (A to C)
        pts2 <- matrix(rep(C, length(indented)), ncol = 2, byrow = TRUE) + indented * (matrix(rep(A, length(indented)), ncol = 2, byrow = TRUE) - matrix(rep(C, length(indented)), ncol = 2, byrow = TRUE))
        ## axis 3 (B to C)
        pts3 <- matrix(rep(B, length(indented)), ncol = 2, byrow = TRUE) + indented * (matrix(rep(C, length(indented)), ncol = 2, byrow = TRUE) - matrix(rep(B, length(indented)), ncol = 2, byrow = TRUE)) 
        
        object@s.misc$lgrid <- list(pts1 = pts1, pts2 = pts2, pts3 = pts3, posgrid = vdivision)
        assign(name_obj, object, envir = parent.frame())
    })


setMethod(
    f = "panelbase",
    signature = "ADEg.Tr",
    definition = function(object, x, y) {
        callNextMethod()
        ## draw triangle (A -> B , B -> C,  C -> A)
        ## small triangle: points distribution
        
        ## triangle vertices
        dfcorner <- rbind(object@s.misc$cornerp$A, object@s.misc$cornerp$B, object@s.misc$cornerp$C, object@s.misc$cornerp$A)
        
        panel.polygon(dfcorner, col = object@adeg.par$pbackground$col, border = if(object@adeg.par$pbackground$box) col = "#000000" else "transparent") ## not really useful (only  for arguments consistency)

        ## size of the grid
        nn <- sapply(object@s.misc$lgrid, nrow)[-4]
        ## draw grid
        if(object@adeg.par$pgrid$draw)
            panel.segments(x0 = c(rep(object@s.misc$lgrid[[1L]][, 1], 2), object@s.misc$lgrid[[2L]][, 1]), 
                           x1 = c(rev(object@s.misc$lgrid[[2L]][, 1]), rep(rev(object@s.misc$lgrid[[3L]][, 1]), 2)), 
                           y0 = c(rep(object@s.misc$lgrid[[1L]][, 2], 2),object@s.misc$lgrid[[2L]][, 2]), 
                           y1 = c(rev(object@s.misc$lgrid[[2L]][, 2]), rep(rev(object@s.misc$lgrid[[3L]][, 2]), 2)),
                           lwd = object@adeg.par$pgrid$lwd,
                           col = object@adeg.par$pgrid$col,
                           lty = object@adeg.par$pgrid$lty)

        ## draw axes
        axis.text2 <- list()
        axis.text <- trellis.par.get("axis.text")
        axis.text2[c("cex", "col")] <- object@adeg.par$pgrid$text[c("cex", "col")]
        division <- object@s.misc$lgrid$posgrid[-c(1, length(object@s.misc$lgrid$posgrid))]
        pos <- c(1, 3, 3)
        srt <- c(0, 60, -60)

        ## get axes names
        if(object@data$storeData)
            axisN <- colnames(object@data$dfxyz)[c(2, 1, 3)]
        else
            axisN <- colnames(eval(object@data$dfxyz, envir = sys.frame(object@data$frame)))[c(2, 1, 3)]
        
        lab <- apply(object@s.misc$lgrid$posgrid, 2, as.character)
        labels <- lab[-c(1, nrow(lab)), ] ## without corner
        
        ## final limits for axes
        lcorners <- lab[c(1, nrow(lab)), ] ## corner lab (limits)
        orderCplot <- dfcorner[c(3, 1, 1, 2, 2, 3), ] ## ordre dessin label, selon row de dfcorner, a reprendre
        posCplot <- rep(c(2, 1, 4), each = 2)
        order_lab <- c(2, 1, 3)
        for(i in 1:3) { ## for the three axis
            ## ticks
            if(object@adeg.par$paxes$draw)
                do.call("panel.text", c(list(labels = labels[, order_lab[i]], x = object@s.misc$lgrid[[i]][, 1], y = object@s.misc$lgrid[[i]][, 2], pos = pos[i], srt = srt[i]), axis.text2))
            ptlab <- object@s.misc$lgrid[[i]][1, ] + (object@s.misc$lgrid[[i]][nn[i], ] - object@s.misc$lgrid[[i]][1, ]) / 2
            
            ## axis names
            do.call("panel.text", args = c(list(labels = axisN[i], x = ptlab[1], y = ptlab[2], srt = srt[i], pos = pos[i]), axis.text))
        }
        do.call("panel.text", c(list(x = orderCplot[, 1], y = orderCplot[, 2], lab = lcorners, pos = posCplot), axis.text2))
    })


setMethod(
    f = "gettrellis",
    signature = "ADEg.Tr",
    definition = function(object) {
        tmp_trellis <- do.call(what = object@lattice.call$graphictype, args = c(formula(1 ~ 1), object@lattice.call$arguments, environment()))
        return(tmp_trellis)
    })


setMethod(
    f = "setlatticecall", 
    signature = "ADEg.Tr",
    definition = function(object) {
        name_obj <- deparse(substitute(object))

        ## background and box
        ## object@trellis.par$panel.background$col <- object@adeg.par$pbackground$col
        if(!object@adeg.par$pbackground$box)
            object@trellis.par$axis.line$col <- "transparent"

        arguments = list(
            par.settings = object@trellis.par,
            scales = if(!is.null(object@g.args$scales)) object@g.args$scales else list(draw = FALSE),
            key = createkey(object),
            aspect = object@adeg.par$paxes$aspectratio,
            panel = function(...) {
                panelbase(object, ...)
                panel(object, ...)
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
