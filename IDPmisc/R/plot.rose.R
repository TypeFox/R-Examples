## plot.rose.R


plot.rose <-
    function (x,
              transf = function(x) sqrt(x),
              subset.col = NULL,
              warn = TRUE,
              general = general.control(),
              grid = grid.control(),
              title = title.control(),
              key = key.control())
    ## Author: Rene Locher
    ## Version: 2009-03-13
    ##
    ## NA are plotted as 0

{
    if (!is.unit(general$mar))
        general$mar <- unit.c(unit(general$mar*grid$cyclVar$cex,
                                   units="char"))
    if (!is.unit(key$between)) key$between <-
        unit(key$between*grid$cyclVar$cex, units="char")

    ## convert distances in viewport independent units
    ## -> distroyes scalability of plots!
    pushViewport(viewport(gp=gpar(cex=general$cex)),recording=FALSE)
    general$mar <- convertWidth(general$mar,"mm")
    key$between <- convertWidth(key$between,"mm")
    popViewport()

    if (!is.null(key$x)&&!is.unit(key$x))
        key$x <- unit(key$x,units="snpc")

    ## checking if length of cyclVar$lab is compatible with ray$n
    if ((grid$ray$n %% length(grid$cyclVar$lab)))
        stop("'ray.n' is not a multiple of length of 'cyclVar.lab'!\n")

    rho <- x@rho
    if (!is.null(subset.col)) { # chose columns to be plotted
        if (is.numeric(subset.col)) {
            if (any(cc <- !is.element(subset.col,1:ncol(rho))))
                stop(paste("", paste(subset.col[cc],collapse=", "),
                  "in 'subset.col' define(s) no valid column of 'x'!"))
        } else
        if (any(cc <- !is.element(subset.col,colnames(rho))))
            stop(paste("", paste(subset.col[cc],collapse=", "),
                   "in 'subset.col' define(s) no valid column of 'x'!"))
        rho <- rho[,subset.col,drop=FALSE]
    }

    ## check if transformation function is compatible with data
    if (any(!is.finite(transf(as.vector(rho))) &
            !is.na(as.vector(rho))))
      stop("Transformation function 'transf' is incompatible with data. You possibly have tried to apply the square root transformation to negativ data.\n")

    ## make sure that stacked and unstacked case have
    ## the same order in the legend
    if (!general$stacked) general$rev.col <- !general$rev.col
    if (general$rev.col) rho <- rho[,ncol(rho):1,drop=FALSE]

    ## NAs are interpreted as 0 (when counts or percentages)
    ## in other cases this is arbitrary
    i.na <- is.na(rho)
    if (sum(i.na)) {
        rho[i.na] <- 0
        if (warn) warning("NAs encountered in rho! NAs are set to 0. Be careful with interpretation!\n",
                          call. = FALSE)
    }

    if (general$stacked) {
        if (min(rho,na.rm=TRUE)<0)
            stop("Stacked roses make sense only for positive variables like counts,  proportions and concentrations!")
        if (ncol(rho)>1)
            rho <- t(apply(rho, MARGIN=1, cumsum))[,ncol(rho):1]
    } ## stacked

    ## calculating the labels for the (main) circles
    if (is.null(grid$circ$r)) {
        if (is.null(grid$ray$lim)) {
            grid$circ$value <- pretty(c(0,rho), n = grid$circ$n)
        } else
        grid$circ$value <- pretty(c(grid$ray$lim, n = grid$circ$n))
    } else grid$circ$value <- grid$circ$r

    nc <- ncol(rho)

    ## choosing adequate colors for displaying the data
    if(general$stacked){ ## well blended colors
        if (is.null(general$col)) general$col <-
            IDPcolorRamp(nc,
                         colInt = data.frame(
                         h = c(0.6, 0.55, 0.45, 0.25),
                         s = c(0.5, 0.55, 0.55, 0.55),
                         v = c(0.92, 0.92, 0.92,0.92)),
                         fr     = c(0.4,0.3))
    } else { ## distinct colors
        if (is.null(general$col))
            general$col <-
                c("#324B80", "#198019", "#DB1919", "#ED9900", "#6E286E")
##                 c("blue","green4" ,"red", "darkorchid4", "black",
##                   "deepskyblue","green","orange", "violetred",
##                   "grey50", "saddlebrown")
    }

    ## constructing a col vector of correct length
    ll <- length(general$col)%/%nc + 1
    if (length(general$col)>nc) general$col <- rep(general$col,ll)[1:nc]

    ## constructing a lty vector of correct length
    ll <- length(general$lty)%/%nc + 1
    if (length(general$lty)>nc) general$lty <- rep(general$lty,ll)[1:nc]

    ## define ray.lim if NULL or
    ## trim circ$value if ray.lim defines smaller range
    if (is.null(grid$ray$lim))
        grid$ray$lim <- range(grid$circ$value) else
    grid$circ$value <-
        grid$circ$value[grid$circ$value >= grid$ray$lim[1] &
                        grid$circ$value<=grid$ray$lim[2]]
    if (grid$ray$lim[1]!=0 & warn)
        warning("Be careful! The center of the rose is not 0, which might be misleading to the reader.\n")

    ## calculating circle radius in transformated native coordinates:
    ## center of circles is trans(ray.lim[1])
    ## untransformed radii are retained in circ$value
    grid$circ$r <- transf(grid$circ$value) - transf(grid$ray$lim[1])

    ## adjust circ$n to the actual number of circles
    grid$circ$n <- length(grid$circ$r)

    ## subcircles are plotted when either circ$sub$n or circ$sub$rle
    ## are defined
    grid$circ$sub$plot <- (!is.null(grid$circ$sub$n) |
                           !is.null(grid$circ$sub$r))

    if (grid$circ$sub$plot){
        if (length(grid$circ$value)>1) {
            if (is.null(grid$circ$sub$r))
                grid$circ$sub$r <-
        transf(seq(grid$circ$value[1],
                   by = (diff(grid$circ$value[1:2]) / grid$circ$sub$n),
                   to = grid$circ$value[length(grid$circ$value)])) -
                        transf(grid$ray$lim[1]) else
            grid$circ$sub$r <-
                transf(grid$circ$sub$r) - transf(grid$ray$lim[1])
        } else {
            grid$circ$sub$plot <- FALSE
            if (warn) warning("Definitions of 'circ.n', 'circ.r' incompatible with definitions for subcircles. No subcircles are plotted.\n")
        }
    }

    if (key$plot&&nc>1) {## plot legend
        if (is.null(key$lab)) key$lab <- colnames(rho)

        if (general$stacked)
            key.grob <- draw.leg(key = list(rectangles = list(
                                            col=general$col,
                                            size = 2,
                                            lwd = 0.5),
                                 text = list(key$lab),
                                 cex = general$cex,
                                 between = 1,
                                 between.rows = 0.5,
                                 between.title = 0.5*grid$cyclVar$cex,
                                 title = key$title,
                                 adj.title = 0,
                                 cex.title = grid$cyclVar$cex,
                                 transparent = TRUE))
        else
            key.grob <- draw.leg(key = list(lines = list(
                                            col = general$col,
                                            lwd = general$lwd,
                                            lty = general$lty,
                                            size = 3),
                                 text = list(key$lab),
                                 between = 1,
                                 between.rows = 0.5,
                                 between.title = 0.5*grid$cyclVar$cex,
                                 cex = general$cex,
                                 title = key$title,
                                 adj.title = 0,
                                 cex.title = grid$cyclVar$cex,
                                 transparent = TRUE))

        ## convert distances in viewport independent units
        pushViewport(viewport(gp=gpar(cex=general$cex)),recording=FALSE)
        keyWidth <- grobWidth(key.grob)
        if (!is.null(key$title))
            keyWidth <- convertWidth(max(grobWidth(key.grob),
                                         unit(grid$cyclVar$cex,
                                              "strwidth",key$title)),"mm") else
        keyWidth <- convertWidth(grobWidth(key.grob),"mm")

        popViewport()
    } else
    keyWidth <- unit(0,"mm")

    ##------------

    if (general$rose$auto) {
        ## calculate proper rose$rad, when it is not defined explicitely
        delta <- unit(max(sapply(grid$cyclVar$lab,nchar)),"lines")
        general$rose$rad <- unit(0.5,"snpc") - 0.5*keyWidth -
            0.5*key$between - delta
        rose.x <-  rose.y <- general$rose$rad+delta
    } else delta <- unit(max(sapply(grid$cyclVar$lab,nchar)),"lines")

    if (is.null(general$rose$x)) rose.x <- general$rose$rad+delta else
    rose.x <- general$rose$x

    if (is.null(general$rose$y)) rose.y <- general$rose$rad+delta else
    rose.y <- general$rose$y

    vp.rose <- viewport(name="vp.rose",
                        x = rose.x,
                        y = rose.y,
                        width = 2*general$rose$rad,
                        height = 2*general$rose$rad,
                        xscale = c(transf(grid$ray$lim)[1]-
                        transf(grid$ray$lim)[2],
                        transf(grid$ray$lim)[2]-
                        transf(grid$ray$lim)[1]),
                        yscale = c(transf(grid$ray$lim)[1]-
                        transf(grid$ray$lim)[2],
                        transf(grid$ray$lim)[2]-
                        transf(grid$ray$lim)[1]),
                        just = c("center","center"),
                        gp = gpar(cex=general$cex),
                        clip = "off")

    gd <- griddat(rho = rho,
                  cyclVar = x@cyclVar,
                  circle = x@circle,
                  vp = vp.rose,
                  grid = grid,
                  title = title)

    if (general$rose$auto) {
        ## readjustment of rose$rad
        if (convertWidth(key$between,"mm",valueOnly=TRUE)>0)
            delta <- 0.5*key$between else delta <- 0*key$between

        rad1 <- convertWidth(0.5*unit(1,"npc") - 0.5*keyWidth -
                             0.5*sum(gd$labSpace[c(2,4)]) -
                             0.5*sum(general$mar[c(2,4)]) - delta, "mm")
        rad2 <- convertHeight(0.5*unit(1,"npc") -
                              0.5*sum(gd$labSpace[c(1,3)]) -
                              0.5*sum(general$mar[c(1,3)]) -
                              delta, "mm")

        general$rose$rad <- min(unit.c(rad1,rad2))
    }

    if (is.null(general$rose$x))
        rose.x <- general$rose$rad+gd$labSpace[2] + general$mar[2] else
    rose.x <- general$rose$x

    if (is.null(general$rose$y))
        rose.y <- general$rose$rad+gd$labSpace[1] + general$mar[1] else
    rose.y <- general$rose$y

    vp.rose <- viewport(name="vp.rose",
                        x = rose.x,
                        y = rose.y,
                        width = 2*general$rose$rad,
                        height = 2*general$rose$rad,
                        xscale = c(transf(grid$ray$lim)[1]-
                        transf(grid$ray$lim)[2],
                        transf(grid$ray$lim)[2]-
                        transf(grid$ray$lim)[1]),
                        yscale = c(transf(grid$ray$lim)[1]-
                        transf(grid$ray$lim)[2],
                        transf(grid$ray$lim)[2]-
                        transf(grid$ray$lim)[1]),
                        just = c("center","center"),
                        gp = gpar(cex=general$cex),
                        clip = "off")
    pushViewport(vp.rose)
    grid.draw(rose.grob(rho = rho,
                        cyclVar = x@cyclVar,
                        circle = x@circle,
                        transf = transf,
                        general = general,
                        grid = grid,
                        title = title,
                        gdat = gd))
    upViewport()

    if (key$plot&&nc>1) {## plot legend
        vp.key <-  viewport(x = if (is.null(key$x))
                            2*general$rose$rad + general$mar[2] +
                            sum(gd$labSpace[c(2,4)]) +
                            key$between else key$x,
                            y = general$mar[1],
                            default.units = "npc",
                            width = grobWidth(key.grob),
                            height = grobHeight(key.grob),
                            just = c("left","bottom"),
                            name = "vp.key")
        pushViewport(vp.key)
        grid.draw(key.grob)
        upViewport()
    }
} ## plot.rose
