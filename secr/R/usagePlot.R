usagePlot <- function(object, add = FALSE, occasion = NULL, col =
                       'black', fill = FALSE, scale = 2, metres = TRUE, rad = 5, ...) {

    if (ms(object))  {
        lapply(object, usagePlot, add = add, occasion = occasion, col = col,
               fill = fill, scale = scale, metres = metres, rad = rad, ...)
    }
    else {
        if (is.null(usage(object)))
            stop ("object does not have usage attribute")
        if (!add)
            plot(object, ...)

        nocc <- ncol(usage(object))
        if (is.null(occasion)) {
            dx  <- rep((cos((1:nocc) * 2 * pi / nocc) * rad), each=nrow(object))
            dy  <- rep((sin((1:nocc) * 2 * pi / nocc) * rad), each=nrow(object))
            xy <- cbind(rep(object$x, nocc) + dx, rep(object$y, nocc)-dy)
            radius <- as.numeric(sqrt(usage(object)) * scale)
        }
        else {
            xy <- object
            if (tolower(occasion[1]) == "all")
                occasion <- 1:nocc
            usge <- apply(usage(object)[,occasion, drop = FALSE], 1, sum)
            radius <- sqrt(usge) * scale
        }

        # metres: use symbols with inches = FALSE
        if (metres) {
            fg <- col
            if (fill)
                bg <- fg
            else
                bg <- NULL
            symbols(xy, circles = radius, inches = FALSE, fg = fg, bg = bg, add = T)
        }
        else {
            pch <- ifelse (fill, 16, 1)
            points(xy, cex = radius, pch = pch, col = col)
        }
    }
}
