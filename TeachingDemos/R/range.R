`%<%` <- function(x,y) {
	xx <- attr(x,'orig.y')
	yy <- attr(y,'orig.x')

	if(is.null(xx)) {
		xx <- x
		x <- rep(TRUE, length(x))
	}
	if(is.null(yy)) {
		yy <- y
		y <- rep(TRUE, length(y))
	}

	out <- x & y & (xx < yy)
	attr(out, 'orig.x') <- xx
	attr(out, 'orig.y') <- yy
	
	out
}

`%<=%` <- function(x,y) {
	xx <- attr(x,'orig.y')
	yy <- attr(y,'orig.x')

	if(is.null(xx)) {
		xx <- x
		x <- rep(TRUE, length(x))
	}
	if(is.null(yy)) {
		yy <- y
		y <- rep(TRUE, length(y))
	}

	out <- x & y & (xx <= yy)
	attr(out, 'orig.x') <- xx
	attr(out, 'orig.y') <- yy
	
	out
}




x <- -3:3

 -2 %<% x %<% 2 
c( -2 %<% x %<% 2 )
x[ -2 %<% x %<% 2 ]
x[ -2 %<=% x %<=% 2 ]


x <- rnorm(100)
y <- rnorm(100)

x[ -1 %<% x %<% 1 ]
range( x[ -1 %<% x %<% 1 ] )


cbind(x,y)[ -1 %<% x %<% y %<% 1, ]
cbind(x,y)[ (-1 %<% x) %<% (y %<% 1), ]
cbind(x,y)[ ((-1 %<% x) %<% y) %<% 1, ]
cbind(x,y)[ -1 %<% (x %<% (y %<% 1)), ]
cbind(x,y)[ -1 %<% (x %<% y) %<% 1, ] # oops













