## latticist: a Lattice-based exploratory visualisation GUI
##
## Copyright (c) 2008 Felix Andrews <felix@nfrac.org>
## GPL version 2 or newer

panel.usertext <-
    function(x, y = NULL, labels = seq_along(x), col = user.text$col,
             alpha = user.text$alpha, cex = user.text$cex, srt = 0, lineheight = user.text$lineheight,
             font = user.text$font, fontfamily = user.text$fontfamily, fontface = user.text$fontface,
             adj = c(0.5, 0.5), pos = NULL, offset = 0.5, ...)
{
    user.text <- current.user.text()
    panel.text(x, y, labels, col = col, alpha = alpha, cex = cex, srt = srt,
               lineheight = lineheight, font = font, fontfamily = fontfamily,
               fontface = fontface, adj = adj, pos = pos, offset = offset, ...)
}

current.user.text <- function() {
    user.text <- trellis.par.get("user.text")
    if (is.null(eval(user.text))) {
        user.text <- trellis.par.get("add.text")
    }
    user.text
}

current.brush.symbol <- function() {
    brush.symbol <- trellis.par.get("brush.symbol")
    if (is.null(eval(brush.symbol))) {
        ## defaults:
        brush.symbol <-
            list(pch = 21, col = "black", fill = "yellow",
                 alpha = 1, cex = 0.8, font = 1)
        plot.symbol <- trellis.par.get("plot.symbol")
        ## take cex from plot.symbol
        brush.symbol$cex <- plot.symbol$cex
        ## use filled equivalent to current plot symbol
        ppch <- as.character(plot.symbol$pch)
        pch <- 21 ## circle, default
        if (ppch %in% c("0", "7", "12", "15", "22")) {
            pch <- 22 ## square
        } else if (ppch %in% c("5", "9", "18", "23")) {
            pch <- 23 ## diamond
        } else if (ppch %in% c("2", "17", "24")) {
            pch <- 24 ## up triangle
        } else if (ppch %in% c("6", "25")) {
            pch <- 25 ## down triangle
        }
        brush.symbol$pch <- pch
    }
    brush.symbol
}

current.brush.line <- function() {
    brush.line <- trellis.par.get("brush.line")
    if (is.null(eval(brush.line))) {
        ## defaults:
        brush.line <- list(col = "red", alpha = 1,
                           lwd = 2, lty = 1)
    }
    brush.line
}

panel.brushpoints <-
    function(x, y = NULL, col = brush.symbol$col, pch = brush.symbol$pch,
             alpha = brush.symbol$alpha, fill = brush.symbol$fill, cex = brush.symbol$cex, ...)
{
    brush.symbol <- current.brush.symbol()
    panel.points(x, y, col = col, pch = pch, alpha = alpha,
                 fill = fill, cex = cex, ...)
}

panel.brushlines <-
    function(x, y = NULL, type = "l", col = brush.line$col,
             alpha = brush.line$alpha, lty = brush.line$lty,
             lwd = brush.line$lwd, ...)
{
    brush.line <- current.brush.line()
    panel.lines(x, y, type = type, col = col, alpha = alpha,
                lty = lty, lwd = lwd, ...)
}
