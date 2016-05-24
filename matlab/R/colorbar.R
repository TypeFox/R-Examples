###
### $Id: colorbar.R 48 2014-02-05 20:50:54Z plroebuck $
###
### Displays colorbar showing the color scale.
###

SIDE.BELOW <- as.integer(1)
SIDE.LEFT  <- as.integer(2)
SIDE.ABOVE <- as.integer(3)
SIDE.RIGHT <- as.integer(4)

LAS.PARALLEL      <- as.integer(0)
LAS.HORIZONTAL    <- as.integer(1)
LAS.PERPENDICULAR <- as.integer(2)
LAS.VERTICAL      <- as.integer(3)


##-----------------------------------------------------------------------------
colorbar <- function(C, location = c("EastOutside",
                                     "WestOutside",
                                     "NorthOutside",
                                     "SouthOutside"), ...) {
    if (!is.numeric(C)) {
        stop(sprintf("argument %s must be numeric", sQuote("C")))
    }

    if (!is.character(location)) {
        stop(sprintf("argument %s must be character", sQuote("location")))
    }
    location <- match.arg(location)

    range.C <- range(C)
    min.lim.C <- round(range.C[1])
    max.lim.C <- round(range.C[2])
    seq.lim.C <- seq(min.lim.C, max.lim.C)

    colorbar.EO <- function(...) {
        ## Swap number of lines of margin on left and right sides
        ## to give enough room to draw the tickmark labels on right
        saved.par <- par(mar = c(5, 2, 4, 4) + 0.1)
        on.exit(par(saved.par))

        graphics::image(x    = 1,
                        y    = seq.lim.C,
                        z    = t(matrix(seq.lim.C, ncol = 1)),
                        ann  = FALSE,
                        xaxt = "n",
                        yaxp = c(min.lim.C, max.lim.C, 5),
                        axes = FALSE,
                        ...)
        axis(SIDE.RIGHT, axTicks(SIDE.RIGHT), las = LAS.HORIZONTAL)
        box()
    }

    colorbar.WO <- function(...) {
        graphics::image(x    = 1,
                        y    = seq.lim.C,
                        z    = t(matrix(seq.lim.C, ncol = 1)),
                        ann  = FALSE,
                        xaxt = "n",
                        yaxp = c(min.lim.C, max.lim.C, 5),
                        las  = LAS.HORIZONTAL,
                        ...)
    }

    colorbar.NO <- function(...) {
        graphics::image(x    = seq.lim.C,
                        y    = 1,
                        z    = t(matrix(seq.lim.C, nrow = 1)),
                        ann  = FALSE,
                        xaxp = c(min.lim.C, max.lim.C, 5),
                        yaxt = "n",
                        axes = FALSE,
                        ...)
        axis(SIDE.ABOVE, axTicks(SIDE.ABOVE))
        box()
    }

    colorbar.SO <- function(...) {
        graphics::image(x    = seq.lim.C,
                        y    = 1,
                        z    = t(matrix(seq.lim.C, nrow = 1)),
                        ann  = FALSE,
                        xaxp = c(min.lim.C, max.lim.C, 5),
                        yaxt = "n",
                        ...)
    }

    switch(EXPR = location,
           EastOutside  = colorbar.EO(...),
           WestOutside  = colorbar.WO(...),
           SouthOutside = colorbar.SO(...),
           NorthOutside = colorbar.NO(...))

    invisible(NULL)
}

