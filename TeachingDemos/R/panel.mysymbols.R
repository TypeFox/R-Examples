panel.my.symbols <- function(x, y, symb, inches=1, polygon = FALSE,
                             ..., symb.plots=FALSE, subscripts, MoreArgs ) {
    if(symb.plots) {
        stop('self plotting symbols (symb.plots=TRUE) is not implemented yet')
    }

    dots <- list(...)
    tmp <- sapply(dots, is.null)
    dots[tmp] <- NULL
    if ( 'type' %in% names(dots) ) dots$type <- 'l'
    tmp.xlen <- length(x)

    if( (length(inches) != 1) && (length(inches) != tmp.xlen) ) {
        inches <- rep(inches[subscripts], length.out=tmp.xlen)
    }

    dots <- lapply(dots, function(x) {
        if( (length(x) != 1) && (length(x) != tmp.xlen) )  {
            x <- rep(x[subscripts], length.out=tmp.xlen)
        }
        x } )

    plotfun <- if( is.function(symb) ) {
        function(x,y,inches,polygon,symb, ...) {
            dots1 <- list(...)
            sargs <- setdiff(names(formals(symb)),'...')
            dots2 <- dots1[sargs]
            dots1[sargs] <- NULL
            symb2 <- xy.coords(do.call(symb,dots2))
            xx <- grid::convertWidth( grid::unit(symb2$x*inches/2, 'inches'),
                                 'native', TRUE )
            yy <- grid::convertHeight( grid::unit(symb2$y*inches/2, 'inches'),
                                 'native', TRUE )
            dots1$x <- x+xx
            dots1$y <- y+yy
            if(polygon) {
                do.call(lattice::lpolygon, dots1)
            } else {
                do.call(lattice::llines, dots1)
            }
        }
    } else {
        function(x,y,inches,polygon,symb, ...) {
            dots <- list(...)
            symb2 <- xy.coords(symb)
            xx <- grid::convertWidth( grid::unit(symb2$x*inches/2, 'inches'),
                                 'native', TRUE )
            yy <- grid::convertHeight( grid::unit(symb2$y*inches/2, 'inches'),
                                 'native', TRUE )
            dots$x <- x+xx
            dots$y <- y+yy
            if(polygon) {
                do.call(lattice::lpolygon, dots)
            } else {
                do.call(lattice::llines, dots)
            }
        }
    }

    funargs <- c(list(x=x, y=y, inches=inches, polygon=polygon),
                 dots)
    funargs$FUN <- plotfun
    if(missing(MoreArgs)) {
        funargs$MoreArgs <- list(symb=symb)
    } else {
        funargs$MoreArgs <- c(MoreArgs, list(symb=symb))
    }

    do.call(mapply, funargs)

    invisible(NULL)

}





### original code
if(FALSE) {

my.df <- data.frame( x=runif(10), y=runif(10) )

xyplot(y~x, my.df, panel=function(x,y,...) {
	xx <- grid::convertX( grid::unit(ms.male[,1]/5, 'inches'), 'native', TRUE )
	yy <- grid::convertY( grid::unit(ms.male[,2]/5, 'inches'), 'native', TRUE )

	xx <- c(xx,NA); yy <- c(yy, NA)

	llines( outer(xx, x, '+'), outer(yy, y, '+') )
	}

)

}

# convert and unit from grid package
