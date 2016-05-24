addpoly.rma <-
function (x, row = -2, level = x$level, digits = 2, annotate = TRUE, 
    mlab, transf, atransf, targs, efac = 1, col, border, cex, 
    ...) 
{
    if (!is.element("rma", class(x))) 
        stop("Argument 'x' must be an object of class \"rma\".")
    if (!x$int.only) 
        stop("Fitted model should not contain moderators.")
    if (missing(mlab)) 
        mlab <- NULL
    if (missing(transf)) 
        transf <- FALSE
    if (missing(atransf)) 
        atransf <- FALSE
    if (missing(targs)) 
        targs <- NULL
    if (missing(cex)) 
        cex <- NULL
    if (missing(col)) 
        col <- "black"
    if (missing(border)) 
        border <- "black"
    if (is.null(mlab)) 
        mlab <- ifelse((x$method == "FE"), "FE Model", "RE Model")
    addpoly(x$b, ci.lb = x$ci.lb, ci.ub = x$ci.ub, rows = row, 
        level = level, digits = digits, annotate = annotate, 
        mlab = mlab, transf = transf, atransf = atransf, targs = targs, 
        efac = efac, col = col, border = border, cex = cex, ...)
}
