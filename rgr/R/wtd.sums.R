wtd.sums <-
function (xx, ri, xloc = NULL, xspread = NULL) 
{
    if (!is.matrix(xx)) 
        stop(deparse(substitute(xx)), " is not a Matrix")
    temp.x <- remove.na(xx)
    x <- temp.x$x
    ncolx <- temp.x$m
    nlnri <- length(ri)
    if (ncolx != nlnri) 
        stop("\n  Number of variables and importances do not match")
    w <- ri/sum(abs(ri))
    a <- w/sqrt(sum(w * w))
    if (is.null(xloc)) 
        xloc <- apply(x, 2, median)
    else {
        if (length(xloc) != nlnri) 
            stop("\n  Numbers of variables and locations do not match")
    }
    if (is.null(xspread)) 
        xspread <- apply(x, 2, mad)
    else {
        if (length(xspread) != nlnri) 
            stop("\n  Numbers of variables and spreads do not match")
    }
    xcentred <- sweep(x, 2, xloc, "-")
    z <- sweep(xcentred, 2, xspread, "/")
    ws <- z %*% a
    invisible(list(input = deparse(substitute(xx)), xloc = xloc, 
        xspread = xspread, ri = ri, w = w, a = a, ws = ws))
}
