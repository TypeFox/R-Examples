HWidentify <- function(x,y,label=seq_along(x), lab.col='darkgreen',
                           pt.col='red', adj=c(0,0), clean=TRUE,
                           xlab=deparse(substitute(x)),
			   ylab=deparse(substitute(y)),  ...) {

	plot(x,y,xlab=xlab, ylab=ylab,...)

	dx <- grconvertX(x,to='ndc')
	dy <- grconvertY(y,to='ndc')

	mm <- function(buttons, xx, yy) {
		d <- (xx-dx)^2 + (yy-dy)^2
		if ( all( d > .01 ) ){
			plot(x,y,xlab=xlab,ylab=ylab,...)
			return()
		}
		w <- which.min(d)
		plot(x,y,xlab=xlab,ylab=ylab,...)
		points(x[w],y[w], cex=2, col=pt.col)
		text(grconvertX(xx,from='ndc'),grconvertY(yy,from='ndc'),
			label[w], col=lab.col, adj=adj)
		return()
	}

	md <- function(buttons, xx, yy) {
		if (any(buttons=='2')) return(1)
		return()
	}

	getGraphicsEvent('Right Click to exit', onMouseMove = mm, onMouseDown=md)
        if(clean) mm( , Inf, Inf )
	invisible()
}

# tmpx <- runif(25)
# tmpy <- rnorm(25)
# HWidentify(tmpx,tmpy,LETTERS[1:25], pch=letters)




HTKidentify <- function(x,y,label=seq_along(x), lab.col='darkgreen',
                            pt.col='red', adj=c(0,0),
                            xlab=deparse(substitute(x)),
			    ylab=deparse(substitute(y)), ...) {

    if( !requireNamespace("tkrplot", quietly=TRUE) ) stop ('tkrplot package is required')

	dx <- numeric(0)
	dy <- numeric(0)

	xx <- yy <- 0

	replot <- function() {
		d <- (xx-dx)^2 + (yy-dy)^2
		if ( all( d > .01 ) ) {
			plot(x,y,xlab=xlab,ylab=ylab,...)
			if( length(dx)==0 ) {
				dx <<- grconvertX(x, to='ndc')
				dy <<- grconvertY(y, to='ndc')
			}
			return()
		}
		w <- which.min(d)
		plot(x,y,xlab=xlab,ylab=ylab,...)
		points(x[w],y[w], cex=2, col=pt.col)
		text(grconvertX(xx,from='ndc'),grconvertY(yy,from='ndc'),
			label[w], col=lab.col, adj=adj)
	}

	tt <- tcltk::tktoplevel()
	img <- tkrplot::tkrplot(tt, replot, hscale=1.5, vscale=1.5)
	tcltk::tkpack(img, side='top')
	iw <- as.numeric(tcltk::tcl("image","width", tcltk::tkcget(img, "-image")))
	ih <- as.numeric(tcltk::tcl("image","height", tcltk::tkcget(img, "-image")))

	cc <- function(x,y) {
		x <- (as.double(x) -1)/iw
		y <- 1-(as.double(y)-1)/ih
		c(x,y)
	}

	mm <- function(x, y) {
		xy <- cc(x,y)
		xx <<- xy[1]
		yy <<- xy[2]
		tkrplot::tkrreplot(img)
	}

	tcltk::tkbind(img, "<Motion>", mm)

	invisible()
}

