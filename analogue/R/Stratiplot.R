## Stratigraphic plots using lattice
`Stratiplot` <- function(x, ...)
    UseMethod("Stratiplot")

`Stratiplot.default` <- function(x, y,
                                 type = "l",
                                 ylab = NULL,
                                 xlab = "",
                                 pages = 1,
                                 rev = TRUE,
                                 ylim,
                                 sort = c("none", "wa", "var"),
                                 svar = NULL,
                                 rev.sort = FALSE,
                                 strip = FALSE,
                                 topPad = 6,
                                 varTypes = "relative",
                                 absoluteSize = 0.5,
                                 zoneNames = NULL,
                                 drawLegend = TRUE,
                                 na.action = "na.omit",
                                 labelAt = NULL,
                                 labelRot = 60,
                                 yticks,
                                 ...) {
    ## inline function for custom axis
    axis.VarLabs <- function(side, ...) {
        if(isTRUE(all.equal(side, "top"))) {
            M <- function(lims) min(lims) + (diff(lims) / 2)
            xlim <- current.panel.limits()$xlim
            if (is.null(labelAt)) {
                at <- M(xlim)
            } else {
                at <- rep(labelAt, 1)
            }
            panel.axis(side = side, outside = TRUE,
                       at = at,
                       tck = 1, line.col = "black",
                       text.col = "black",
                       labels = levels(sx$ind)[which.packet()],
                       rot = labelRot)
        } else {
            axis.default(side = side, ...)
        }
    }

    ## process 'type'
    TYPE <- c("l","h","g","smooth","b","o", "poly","p")
    if(any(!(type %in% TYPE)))
        stop("Invalid 'type' specified")
    n.vars <- NCOL(x)
    ## ylabel
    if(is.null(ylab)) {
        ylab <- deparse(substitute(y))
    }
    ## do we need to sort variables
    if(missing(sort)) {
        sort <- "none"
    }
    sort <- match.arg(sort)
    if((check.var <- (missing(svar) || is.null(svar)))) {
        svar <- y
    }
    ord <- seq_len(n.vars)
    if(sort == "wa") {
        ## sort by
        opt <- optima(x, svar)
        ord <- order(opt)
    } else if(sort == "var") {
        if(check.var) {
            warning("With 'sort = \"var\"', 'svar' not supplied.\nNo sorting applied.")
        } else {
            ord <- order(svar)
        }
    }
    if (rev.sort) {
        ord <- rev(ord)
    }
    ## apply ordering
    x <- x[, ord]

    ## stack the data
    sx <- stack(x)
    sx$ind <- factor(sx$ind, levels = colnames(x)) # add grouping variable

    ## check length of y
    if(!isTRUE(all.equal((leny <- length(y)), (nr <- nrow(sx))))) {
        ## if length(y) == nrow(sx)/n.vars, then expand
        if(isTRUE(all.equal(leny, nr / n.vars))) {
            y <- rep(y, n.vars)
        } else {
            stop("Ambiguous 'length(y)';\nmust be equal to 'nrow(x)' or\n'nrow(x) / number of species'.")
        }
    }

    ## handle NA's, by default NA pass in here fine, but this messes
    ## up line plots or polygons for example
    NAFUN <- match.fun(na.action)
    sx <- NAFUN(sx)
    ## and apply to `y` too
    if(!is.null(NAS <- na.action(sx))) {
        y <- y[-NAS]
    }

    ## plot parameters
    maxy <- max(y, na.rm = TRUE)
    miny <- min(y, na.rm = TRUE)
    ## add padYlim * range as per base graphics - 4% of range
    padY <- 0.04
    if(missing(ylim)) {
        ##diffy <- padY * (maxy - miny)
        diffy <- maxy - miny
        ylim <- c(miny - (padY * diffy), maxy + (padY * diffy))
    } else {
        ## should these be ylim[1] and ylim[2] ???
        minLim <- min(ylim, na.rm = TRUE)
        maxLim <- max(ylim, na.rm = TRUE)
        ## add padY * range as per base graphics
        diffy <- abs(diff(c(minLim, maxLim)))
        ylim <- if(minLim > maxLim) {
                    c(minLim + (padY * diffy), maxLim - (padY * diffy))
                } else {
                    c(minLim - (padY * diffy), maxLim + (padY * diffy))
                }
    }
    ## Reverse the y-axis?
    if (rev) {
        ylim <- rev(ylim)
    }
    ## process the column/variable types
    ## If varTypes of length one replicate it to NCOL(x)
    if(isTRUE(all.equal(length(varTypes), 1L))) {
        varTypes <- rep(varTypes, length = n.vars)
    }
    ## If typeLen != 1 or NCOL shout warning
    if(length(varTypes) != n.vars) {
        warning("Length of 'varTypes' not 1 or equal to number of variables. Recycling or truncating of 'varTypes' as a result.")
    }
    ## Only allow two types of variables: "relative", "absolute"
    if(any(!(varTypes %in% c("relative", "absolute")))) {
        stop("Ambiguous entry in 'varTypes'.\nMust be one of \"relative\", or \"absolute\"")
    }
    ## compute max abundances per relative column, which is used
    ## to scale the panel widths layout.widths parameter)
    max.abun <- sapply(x, function(x) round(max(x, na.rm = TRUE), 1),
                       USE.NAMES = FALSE)
    ## absolute panels should be set to absoluteSize of max.abun
    panelWidths <- max.abun
    ABS <- which(varTypes == "absolute")
    REL <- which(varTypes == "relative")
    if(any(REL)) {
        panelWidths[ABS] <- absoluteSize * max(max.abun[REL], na.rm = TRUE)
    } else {
        panelWidths[ABS] <- absoluteSize
    }
    ## xlim in xyplot call
    xlimits <- lapply(max.abun * 1.05, function(x) c(0, x))
    if(any(ABS)) {
        ## but need any "absolute" panels setting to +/- 0.05(range)
        min.vars <- sapply(x[ABS], min, na.rm = TRUE, USE.NAMES = FALSE)
        max.vars <- sapply(x[ABS], max, na.rm = TRUE, USE.NAMES = FALSE)
        ranges <- (0.04 * (max.vars - min.vars))
        xlimits[ABS] <- as.list(data.frame(t(cbind(min.vars - ranges,
                                                   max.vars + ranges))))
    }
    ## scales in xyplot call
    ## handle custom tick locations for y-axis
    if (missing(yticks)) {
        yticks <- TRUE
    }
    scales <- list(cex = 0.75, tck = 0.75,
                   y = list(axs = "r", limits = ylim, at = yticks),
                   x = list(axs = "r", rot = 45, relation = "free",
                   limits = xlimits))
    par.strip.text <- list(cex = 0.75)
    str.max <- 1

    ## start a new grid page
    grid.newpage()

    if(!isTRUE(strip)) {
        gp <- gpar()
        convWidth <- function(x, gp) {
            convertWidth(grobWidth(textGrob(x, gp = gp)), "lines",
                         valueOnly = TRUE)
        }
        str.max <- max(sapply(levels(sx$ind), convWidth, gp,
                              USE.NAMES = FALSE),
                       na.rm = TRUE)
        str.max <- ceiling(str.max) + topPad
    }
    ## Legend specification for Zones
    dotArgs <- list(...)
    if("zones" %in% names(dotArgs) && drawLegend) {
        Zones <- sort(c(ylim, dotArgs$zones), decreasing = rev)
        Heights <- abs(diff(Zones))
        Ydiff <- abs(diff(ylim))
        Heights <- Heights / Ydiff
        Ylocs <- cumsum(c(0, Heights[-length(Heights)]))
        ZoneUnits <- function(h, units = "npc") {
            unit(h, units = units)
        }
        HeightUnits <- do.call(unit.c, lapply(Heights, ZoneUnits))
        MidUnits <- do.call(unit.c, lapply(Heights/2, ZoneUnits))
        YlocsUnits <- do.call(unit.c, lapply(Ylocs, ZoneUnits))
        ## rectangles to draw
        zoneRects <- rectGrob(y = YlocsUnits, width = unit(1.5, "cm"),
                              height = HeightUnits, vjust = 0)
        if(is.null(zoneNames))
            zoneNames <- paste("Zone", seq_along(c(1, dotArgs$zones)))
        labelText <- rep(zoneNames, length = length(dotArgs$zones) + 1)
        labelYlocs <- YlocsUnits + MidUnits
        zoneLabels <- textGrob(label = labelText, y = labelYlocs,
                               just = rep("center",2),
                               default.units = "npc")
        key.layout <- grid.layout(nrow = 1, ncol = 1,
                                  heights = unit(1, "null"),
                                  widths = unit(1.5, "cm"),
                                  respect = FALSE)
        key.gf <- frameGrob(layout = key.layout)
        key.gf <- placeGrob(key.gf, zoneRects, row = 1, col = 1)
        key.gf <- placeGrob(key.gf, zoneLabels, row = 1, col = 1)
        Legend <- list(right = list(fun = key.gf))
    } else {
        Legend <- NULL
    }
    ## plotting
    plt <- xyplot(y ~ values | ind,
                  data = sx,
                  type = type,
                  ylab = ylab, xlab = xlab,
                  strip.left = FALSE, strip = strip,
                  par.strip.text = par.strip.text,
                  scales = scales,
                  ##xlim = xlimits,
                  ylim = ylim,
                  panel = "panel.Stratiplot",
                  layout = c(n.vars, 1, pages),
                  par.settings = list(layout.widths = list(panel = panelWidths),#max.abun),
                                      layout.heights = list(top.padding = str.max)),
                  axis = if(isTRUE(strip)) {axis.default} else {axis.VarLabs},
                  legend = Legend,
                  ...)
    plot(plt, newpage = FALSE)
    invisible(plt)
}
