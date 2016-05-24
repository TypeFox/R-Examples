siaraddcross <-
function(x = NULL, ex = NULL, y = NULL, ey = NULL,
        clr = "grey50", upch = 21) {
        points(x, y, col = clr, pch = upch)
        if (!is.null(ex)) {
            lines(c(x - ex, x + ex), c(y, y), col = clr)
        }
        if (!is.null(ey)) {
            lines(c(x, x), c(y - ey, y + ey), col = clr)
        }
    }
