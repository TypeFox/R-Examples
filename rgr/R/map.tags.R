map.tags <-
function (xx, yy, tag, xlab = "Easting", ylab = "Northing", taglab = deparse(substitute(tag)), 
    main = "", tol = 0.04, ...) 
{
    frame()
    oldpar <- par()
    on.exit(par(oldpar))
    par(pty = "m")
    temp.x <- remove.na(cbind(xx, yy))
    x <- temp.x$x[1:temp.x$n, 1]
    y <- temp.x$x[1:temp.x$n, 2]
    if (main == "") 
        if (taglab == "") 
            banner <- ""
        else banner <- paste("Map of 'values' for", taglab)
    else banner <- main
    tag[is.na(tag)] <- "+"
    eqscplot(x, y, type = "n", xlab = xlab, ylab = ylab, main = banner, 
        tol = tol, ...)
    text(x, y, tag, ...)
    invisible()
}
