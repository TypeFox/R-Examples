
C_par <- function(x) {
    dev.set(recordDev())
    # Mimic call on off-screen device (so get the right answer when
    # query off-screen device in drawing functions)
    do.call("par", x[-1])
    par <- par()
    dev.set(playDev())
    parnames <- names(x[-1][[1]])
    # Only remake viewports for highest-level change in par()
    if (any(c("oma", "omd", "omi") %in% parnames)) {
        incrementInnerAlpha()
        setUpInner(par)
    } else if (any(c("fig", "fin") %in% parnames)) {
        incrementFigureAlpha()
        setUpFigure(par)
    } else if (any(c("mex", "mai", "mar", "pin", "plt") %in% parnames)) {
        incrementPlotAlpha()
        setUpPlot(par)
    } else if (any(c("usr", "xlog", "ylog") %in% parnames)) {
        # IF we have reset par(usr), we need a new "window" viewport
        incrementWindowAlpha()
        # Align windowPlotAlpha with plotAlpha
        setWindowPlotAlpha(plotAlpha())
        setUpUsr(par$usr)
    }
}

gparParNames <- c("font", "family", "bg", "fg", "col", "lheight",
                  "lend", "ljoin", "lmitre", "ps",
                  "cex", "lex", "lwd", "lty")

gparNameFromParName <- function(x) {
    switch(x,
           font="fontface",
           family="fontfamily",
           bg="fill",
           fg="col",
           lheight="lineheight",
           lend="lineend",
           ljoin="linejoin",
           lmitre="linemitre",
           ps="fontsize",
           x)
}

# 'x' should be a result from calling par() to set new par() values
# (i.e., a list of previous par() values)
gparFromPar <- function(x) {
    gparNames <- sapply(names(x), gparNameFromParName)
    names(x) <- gparNames
    do.call(gpar, x)
}

currentPar <- function(inlinePars) {
    # Drop any inlinePars that are NULL
    # (should never set a par to NULL ?)
    inlinePars <- inlinePars[!sapply(inlinePars, is.null)]
    if (length(inlinePars)) {
        opar <- par(inlinePars)
    }
    cpar <- par()
    # For some weird reason, par()$lty gives different result from par("lty")
    # (and the former can be WRONG or at least an invalid value)
    # Until figure out what is going on, use this workaround
    cpar$lty <- par("lty")
    if (length(inlinePars)) {
        par(opar)
    }
    cpar
}

FixupPch <- function(pch, dflt) {
    if (length(pch) == 0) {
        dflt
    } else {
        pch
    }
}

FixupLty <- function(lty, dflt) {
    if (length(lty) == 0) {
        dflt
    } else {
        lty
    }
}

FixupLwd <- function(lwd, dflt) {
    if (length(lwd) == 0) {
        dflt
    } else {
        ifelse(is.finite(lwd) | lwd >=0, lwd, NA)
    }
}

FixupCol <- function(col, dflt, bg) {
    if (length(col) == 0) {
        dflt
    } else {
        # col=0 means par$bg in 'graphics'
        if (is.numeric(col)) {
            col <- ifelse(col == 0, bg, col)
        }
        col
    }
}

FixupCex <- function(cex, dflt) {
    if (length(cex) == 0) {
        dflt
    } else {
        ifelse(is.finite(cex) & cex > 0, cex, NA)
    }
}

FixupFont <- function(font, dflt) {
    if (length(font) == 0) {
        dflt
    } else {
        if (is.numeric(font)) {
            ifelse(font < 1 | font > 5, NA, font)
        } else {
            font
        }
    }
}

