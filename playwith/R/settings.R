## playwith: interactive plots in R using GTK+
##
## Copyright (c) 2007 Felix Andrews <felix@nfrac.org>
## GPL version 2 or newer

.defaultPlaywithOptions <- function()
    list(
         ## explicit arguments to playwith():
         new = FALSE,
         width = 6,
         height = 6,
         pointsize = 12,
         click.mode = "Zoom",
         time.mode = NA,
         eval.args = NA,
         on.close = NULL,
         ## implicit arguments to playwith():
         page.annotation = FALSE,
         clip.annotations = TRUE,
         label.offset = 0.5,
         arrow = list(length = 0.15, unit = "inches"),
         ## themes:
         themes = alist(
         "WhiteBG" = trellis.par.set(col.whitebg()),
         "Greyscale (for print)" = trellis.par.set(standard.theme("postscript")),
         "DarkBG" = trellis.par.set(standard.theme("X11")),
         "ColorBrewer 1" = trellis.par.set(custom.theme()),
         "ColorBrewer 2" = trellis.par.set(custom.theme.2()),
         "ColorBrewer Black" = trellis.par.set(custom.theme.black())
         ),
         styleShortcuts = alist(
         "Transparent strips" =
         trellis.par.set(list(strip.background = list(col = "transparent"),
                              strip.shingle = list(col = grey(0.1)))),
         "Grey strips" =
         trellis.par.set(list(strip.background = list(col = grey(7:1/8)),
                              strip.shingle = list(col = grey(6:0/8)))),
         "Grey color ramp" =
         trellis.par.set(list(regions = list(col = colorRampPalette(brewer.pal(9, "Greys"))(100)))),
         "Spectral color ramp" =
         trellis.par.set(list(regions = list(col = colorRampPalette(brewer.pal(11, "Spectral"))(100)))),
         "Red vs blue via white" =
         trellis.par.set(list(regions = list(col = colorRampPalette(brewer.pal(11, "RdBu"))(100)))),
         "Dry vs wet (BrBG)" =
         trellis.par.set(list(regions = list(col = colorRampPalette(brewer.pal(11, "BrBG"))(100)))),
         "Yl - Gn - Bu" =
         trellis.par.set(list(regions = list(col = colorRampPalette(brewer.pal(9, "YlGnBu"))(100)))),
         "Yl - Or - Rd" =
         trellis.par.set(list(regions = list(col = colorRampPalette(brewer.pal(9, "YlOrRd"))(100)))),
         "Layout as table" =
         lattice.options(default.args = list(as.table = TRUE))
         ),
         ## global:
         save.as.format = "pdf",
         ui.menus.xml = system.file("etc", "ui.menus.xml", package="playwith"),
         ui.toolbars.xml = system.file("etc", "ui.toolbars.xml", package="playwith"),
         ui.custom.xml = NULL,
         custom.tools = NULL,
         init.actions = NULL,
         preplot.actions = NULL,
         update.actions = NULL,
         custom.toolbar = "BottomToolbar",
         parameters.toolbar = "BottomToolbar",
         toolbar.style = "both",
         show.toolbars = TRUE,
         show.calltoolbar = TRUE,
         show.menubar = TRUE,
         show.statusbar = TRUE,
         show.tooltips = FALSE,
         undo.levels = 4,
         catch.errors = TRUE,
         deparse.options = c()
         )

## code below copied from lattice

playwith.getOption <- function(name)
{
    .PlaywithEnv$options[[name]]
}

playwith.options <- function(...)
{
    ## this would have been really simple if only form allowed were
    ## lattice.options("foo", "bar") and
    ## lattice.options(foo=1, bar=2). But it could also be
    ## lattice.options(foo=1, "bar"), which makes some juggling necessary

    new <- list(...)
    if (is.null(names(new)) && length(new) == 1 && is.list(new[[1]])) new <- new[[1]]
    old <- .PlaywithEnv$options

    ## if no args supplied, returns full options list
    if (length(new) == 0) return(old)

    nm <- names(new)
    if (is.null(nm)) return(old[unlist(new)]) ## typically getting options, not setting
    isNamed <- nm != "" ## typically all named when setting, but could have mix
    if (any(!isNamed)) nm[!isNamed] <- unlist(new[!isNamed])

    ## so now everything has non-"" names, but only the isNamed ones should be set
    ## everything should be returned, however

    retVal <- old[nm]
    names(retVal) <- nm
    nm <- nm[isNamed]

    ## this used to be

    ## modified <- updateList(retVal[nm], new[nm])
    ## .LatticeEnv$lattice.options[names(modified)] <- modified

    ## but then calling lattice.options(foo = NULL) had no effect
    ## because foo would be missing from modified.  So, we now do:

    updateList <- function (x, val) {
        if (is.null(x)) x <- list()
        modifyList(x, val)
    }
    .PlaywithEnv$options <- updateList(old, new[nm])

    ## return changed entries invisibly
    invisible(retVal)
}
