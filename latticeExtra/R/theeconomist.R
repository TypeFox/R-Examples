## Implementation Copyright (c) 2009 Felix Andrews
## based on plot style used in The Economist magazine.

theEconomist.theme <-
   function(win.fontfamily = NULL, #"Gill Sans MT"
            with.bg = FALSE, box = "black", ...)
{
    theme <- list(
         background = list(col = if (with.bg) "#D5E2E9" else "transparent"),
         plot.line = list(col = "#00526D", lwd = 2.5),
         superpose.line = list(col = c("#00526D", "#00A3DB", "#7A2713", "#939598", "#6CCFF6"), lwd = 2.5),
         plot.symbol = list(col = "#00526D", pch = 16),
         superpose.symbol = list(col = c("#00526D", "#00A3DB", "#7A2713", "#939598", "#6CCFF6"), pch = 16),
         plot.polygon = list(col = "#00526D"),
         superpose.polygon = list(col = c("#5F92A8", "#00526D", "#6CCFF6", "#00A3DB", "#A7A9AC")),
         regions = list(col = colorRampPalette(brewer.pal(9, "Blues"))(100)),
         reference.line = list(col = if (with.bg) "white" else "#aaaaaa", lwd = 1.75),
         dot.line = list(col = if (with.bg) "white" else "#aaaaaa", lwd = 1.75),
         add.line = list(col = "#ED1C24", lwd = 1.5),
         axis.line = list(col = box),
         box.3d = list(col = box),
         strip.border = list(col = box),
         strip.background = list(col = if (with.bg) "white" else "#CBDDE6"),
         strip.shingle = list(col = if (with.bg) "#CBDDE6" else "#00A3DB", alpha = 0.5),
         par.main.text = list(font = 1, just = "left", x = grid::unit(5, "mm")),
         par.sub.text = list(font = 1, just = "left", x = grid::unit(5, "mm")),
         axis.text = list(cex = 0.8),
         box.dot = list(col = "#00526D", pch = "|", lwd = 1.75),
         box.rectangle = list(fill = "#00526D", alpha = 0.5, col = "#00526D", lwd = 1.75),
         box.umbrella = list(col = "#00526D", lty = 1, lwd = 1.75)
         )
    if (.Platform$OS.type == "windows" && !is.null(win.fontfamily)) {
        windowsFonts(TheEconomistLike = win.fontfamily)
        theme$grid.pars$fontfamily <- "TheEconomistLike"
    } else {
        ## TODO: how do fonts work on linux etc?
    }
    modifyList(modifyList(standard.theme("pdf"), theme), simpleTheme(...))
}

theEconomist.opts <- function()
{
    list(default.args =
         list(axis = theEconomist.axis,
              xscale.components = xscale.components.subticks,
              between = list(x = 0.8, y = 0.8)),
         axis.padding = list(numeric = 0, factor = 0.6),
         skip.boundary.labels = 0,
         layout.widths =
         list(axis.left = list(x = 0, units = "char"),
              axis.right = list(x = 6, units = "char"))
         )
}

theEconomist.axis <-
    function(side = c("top", "bottom", "left", "right"),
             scales, components, ...,
             labels = c("default", "yes", "no"),
             ticks = c("default", "yes", "no"),
             line.col)
{
    if (scales$draw == FALSE)
        return()
    side <- match.arg(side)
    labels <- match.arg(labels)
    ticks <- match.arg(ticks)
    if (side %in% c("bottom", "top")) {
        if (side == "top")
            ticks <- "no"
        if (scales$relation == "same") {
            scales$alternating <- 1 ## bottom side only
        }
    }
    if (side %in% c("left", "right")) {
        ticks <- "no"
        components[["left"]]$ticks$tck <- 0
        if (scales$relation == "same") {
            scales$alternating <- 2 ## right side only
        } else {
            if (side == "right") {
                labels <- if (scales$draw) "yes" else "no"
                if (!is.list(components$right)) {
                    components$right <- components$left
                }
            }
            if (side == "left") {
                ## check for two different axes on left and right
                if (!is.list(components$right)) {
                    labels <- "no"
                }
            }
        }
    }
    ## use axis.text for ticks because axis.line$col might be transparent
    axis.text <- trellis.par.get("axis.text")
    axis.default(side, scales = scales,
                 components = components, ...,
                 labels = labels, ticks = ticks,
                 line.col = axis.text$col)
    ## now draw grid lines corresponding to horizontal axis ticks.
    ## can only do this with the bottom and right sides;
    ## otherwise the strip viewports are current, not panel.
    if (side %in% c("top", "left"))
        return()
    if (side == "right") {
        comp.list <- components[["right"]]
        if (!is.list(comp.list))
            comp.list <- components[["left"]]
        panel.refline(h = comp.list$ticks$at)
        ## draw axis line along bottom (assuming transparent axis.line)
        lims <- current.panel.limits()
        panel.abline(h = lims$y[1], col = axis.text$col)
    }
}

asTheEconomist <-
    function(x, ...,
             type = "l",
             ylab = expression(NULL),
             xlab = expression(NULL),
             par.settings =
               theEconomist.theme(with.bg = with.bg, box = "transparent"),
             with.bg = FALSE,
             par.strip.text = list(font = 2))
{
    ans <- x
    ## make nice left-aligned title
    title <- ans$main
    if (is.null(title)) title <- ans$ylab
    if (is.null(title)) title <- ans$ylab.default
    ans <- update(ans, main = title,
                  type = type, ylab = ylab, xlab = xlab,
                  par.settings = par.settings,
                  par.strip.text = par.strip.text,
                  between = list(x = 0.8, y = 0.8),
                  scales = list(y = list(axs = "i", alternating = 2)),
                  skip.boundary.labels = 0,
                  lattice.options = list(
                  layout.widths =
                  list(axis.left = list(x = 0, units = "char"),
                       axis.right = list(x = 6, units = "char"))
                  )
                  )
    ## these do not get through update()
    ans$axis <- theEconomist.axis
    ans$xscale.components <- xscale.components.subticks
    ans$call <- match.call()
    ans
}


genGillSans <- function()
{
    ## generate AFM font metrics of GillSans - for ps/pdf
    ## TODO: is there a similar looking free/open font?
    oldwd <- getwd()
    on.exit(setwd(oldwd))
    setwd(system.file("afm", package = "grDevices"))
    fnames <- c("GIL_____", "GILB____", "GILI____", "GILBI___")
    for (fn in fnames) {
        system(sprintf("ttf2afm -o %s.afm %s.TTF", fn, 
                       file.path(Sys.getenv("windir"), "Fonts", fn)))
        system(sprintf("gzip %s.afm", fn))
    }
    TheEconomistLike <- Type1Font("TheEconomistLike",
                                  paste(fnames, ".afm", sep = ""))
    postscriptFonts(TheEconomistLike = TheEconomistLike)
    pdfFonts(TheEconomistLike = TheEconomistLike)
    ps.options(fonts = c("sans", "serif", "mono", "Times", "Helvetica", "Courier", "URWHelvetica", "TheEconomistLike"))
    pdf.options(fonts = ps.options()$fonts)
    ## TODO: problem: after this, plots generated by bitmap() use wrong font
}
