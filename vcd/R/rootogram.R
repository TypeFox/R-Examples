rootogram <- function(x, ...)
{
    UseMethod("rootogram")
}

rootogram.goodfit <- function(x, ...)
{
    rootogram.default(x$observed, x$fitted, names = x$count, df = x$df, ...)
}

rootogram.default <- function(x, fitted, names = NULL, scale = c("sqrt", "raw"),
                              type = c("hanging", "standing", "deviation"),
                              shade = FALSE, legend = TRUE,
                              legend_args = list(x = 0, y = 0.2, height = 0.6), df = NULL,
                              rect_gp = NULL, rect_gp_args = list(), lines_gp = gpar(col = "red", lwd = 2),
                              points_gp = gpar(col = "red"), pch = 19,
                              xlab = NULL, ylab = NULL, ylim = NULL,
                              main = NULL, sub = NULL,
                              margins = unit(0, "lines"), title_margins = NULL,
                              legend_width = NULL,

                              main_gp = gpar(fontsize = 20),
                              sub_gp = gpar(fontsize = 15),

                              name = "rootogram", prefix = "",
                              keep_aspect_ratio = FALSE,
                              newpage = TRUE, pop = TRUE,
                              return_grob = FALSE, ...)
{
    if(is.null(names)) names <- names(x)

    if(is.table(x)) {
        if(length(dim(x)) > 1) stop ("x must be a 1-way table")
        x <- as.vector(x)
    }
    obs <- x
    fit <- fitted
    res <- (obs - fit) / sqrt(fit)
    if(is.null(xlab)) {xlab <-  "Number of Occurrences"}

    if(match.arg(scale) == "sqrt") {
        obs <- sqrt(obs)
        fit <- sqrt(fit)
        if(is.null(ylab)) {ylab <- "sqrt(Frequency)"}
    } else {
        if(is.null(ylab)) {ylab <- "Frequency"}
    }

    ## rect_gp (color, fill, lty, etc.) argument
    if (shade) {
        if (is.null(rect_gp)) rect_gp <- shading_hcl
        if (is.function(rect_gp)) {
            if (is.null(legend) || (is.logical(legend) && legend))
                legend <- legend_resbased
            gpfun <- if (inherits(rect_gp, "grapcon_generator"))
                         do.call("rect_gp", c(list(obs, res, fit, df), rect_gp_args)) else rect_gp
            rect_gp <- gpfun(res)
        } else if (!is.null(legend) && !(is.logical(legend) && !legend))
              stop("rect_gp argument must be a shading function for drawing a legend")
    }

    if (is.null(rect_gp))
        rect_gp <- gpar(fill = "lightgray")

    ## set up page
    if (newpage)
        grid.newpage()
    if (keep_aspect_ratio)
        pushViewport(viewport(width = 1, height = 1, default.units = "snpc"))

    pushViewport(vcdViewport(mar = margins,
                             oma = title_margins,
                             legend = shade && !(is.null(legend) || is.logical(legend) && !legend),
                             main = !is.null(main), sub = !is.null(sub),
                             keep_aspect_ratio = keep_aspect_ratio,
                             legend_width = legend_width,
                             prefix = prefix))

    ## legend
    if (inherits(legend, "grapcon_generator"))
        legend <- do.call("legend", legend_args)
    if (shade && !is.null(legend) && !(is.logical(legend) && !legend)) {
        seekViewport(paste(prefix, "legend", sep = ""))
        legend(res, gpfun, "Pearson\nresiduals:")
    }

    ## titles
    if (!is.null(main)) {
        seekViewport(paste(prefix, "main", sep = ""))
        if (is.logical(main) && main)
            main <- deparse(substitute(x))
        grid.text(main, gp = main_gp)
    }

    if (!is.null(sub)) {
        seekViewport(paste(prefix, "sub", sep = ""))
        if (is.logical(sub) && sub && is.null(main))
            sub <- deparse(substitute(x))
        grid.text(sub, gp = sub_gp)
    }

    seekViewport(paste(prefix, "plot", sep = ""))
    switch(match.arg(type),

           "hanging" = {
               if(is.null(ylim)) {ylim <- range(-0.1 * c(fit-obs,fit),
                                                c(fit-obs,fit)) + c(0, 0.1)}
               dummy <- grid_barplot(obs, names = names, offset = fit - obs, gp = rect_gp,
                                     xlab = xlab, ylab = ylab, ylim = ylim,
                                     name = name, newpage = FALSE, pop = FALSE, return_grob = FALSE, ...)
               downViewport(name)
               grid.lines(x = dummy, y = fit, default.units = "native", gp = lines_gp)
               grid.points(x = dummy, y = fit, default.units = "native", gp = points_gp, pch = pch)
               grid.lines(x = unit(c(0, 1), "npc"), y = unit(0, "native"))
               if(pop) popViewport() else upViewport()
           },

           "standing" = {
               if(is.null(ylim)) {ylim <- range(-0.01 * c(obs,fit), c(obs,fit)) }
               dummy <- grid_barplot(obs, names = names, gp = rect_gp,
                                     xlab = xlab, ylab = ylab, ylim = ylim,
                                     name = name, newpage = FALSE, pop = FALSE, return_grob = FALSE, ...)
               downViewport(name)
               grid.lines(x = dummy, y = fit, default.units = "native", gp = lines_gp)
               grid.points(x = dummy, y = fit, default.units = "native", gp = points_gp, pch = pch)
               if(pop) popViewport() else upViewport()
           },

           "deviation" = {
               if(is.null(ylim)) {ylim <- range(-0.1 * c(fit-obs,fit),
                                                c(fit-obs,fit)) + c(0, 0.1)}
               dummy <- grid_barplot(fit - obs, names = names, gp = rect_gp,
                                     xlab = xlab, ylab = ylab, ylim = ylim,
                                     name = name, newpage = FALSE, pop = FALSE, return_grob = FALSE, ...)
               downViewport(name)
               grid.lines(x = dummy, y = fit, default.units = "native", gp = lines_gp)
               grid.points(x = dummy, y = fit, default.units = "native", gp = points_gp, pch = pch)
               if(pop) popViewport() else upViewport()
           }

           )

    if (return_grob)
        invisible(grid.grab())
    else
        invisible(NULL)
}

grid_barplot <- function(height, width = 0.8, offset = 0,
                         names = NULL, xlim = NULL, ylim = NULL, xlab = "", ylab = "", main = "",
                         gp = gpar(fill = "lightgray"), name = "grid_barplot", newpage = TRUE, pop = FALSE, return_grob = FALSE)
{
    if(is.null(names)) names <- names(height)
    height <- as.vector(height)
    n <- length(height)
    width <- rep(width, length.out = n)
    offset <- rep(offset, length.out = n)

    if(is.null(names)) names <- rep("", n)
    if(is.null(xlim)) xlim <- c(1 - mean(width[c(1, n)]), n + mean(width[c(1, n)]))
    if(is.null(ylim)) ylim <- c(min(offset), max(height + offset))

    if(newpage) grid.newpage()
    pushViewport(plotViewport(xscale = xlim, yscale = ylim, default.units = "native", name = name))
    grid.rect(x = 1:n, y = offset, width = width, height = height,
              just = c("centre", "bottom"), default.units = "native", gp = gp)
    grid.yaxis()
    grid.text(names, x = unit(1:n, "native"), y = unit(rep(-1.5, n), "lines"))
    grid.text(xlab, y = unit(-3.5, "lines"))
    grid.text(ylab, x = unit(-3, "lines"), rot = 90)
    grid.text(main, y = unit(1, "npc") + unit(2, "lines"), gp = gpar(fontface = "bold"))
    if(pop) popViewport() else upViewport()
    if (return_grob)
        invisible(structure(1:n, grob = grid.grab()))
    else
        invisible(1:n)
}
