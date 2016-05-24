updateusr <- function(x1,y1=NULL,x2,y2=NULL) {
    xy1 <- xy.coords(x1,y1)
    xy2 <- if( missing(x2) && missing(y2) ) {
        xy.coords(y1)
    } else {
        xy.coords(x2,y2)
    }

    cur.usr <- par('usr')

    xslope <- diff(xy2$x)/diff(xy1$x)
    yslope <- diff(xy2$y)/diff(xy1$y)

    new.usr.x <- xslope * ( cur.usr[1:2] - xy1$x ) + xy2$x
    new.usr.y <- yslope * ( cur.usr[3:4] - xy1$y ) + xy2$y

    invisible(par(usr=c(new.usr.x, new.usr.y)))
}


# need to add options for dealing with fewer than 2 points, more than 2 points, and NA values.
