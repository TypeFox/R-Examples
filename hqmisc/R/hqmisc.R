# hqmisc 
# version 0.1
# date 2014-03-08

### MISCELLANEOUS FUNCTIONS ###

as.dummies <- function (x) {
	if (!is.factor(x)) { stop("argument must be a factor") }
	result <- matrix( NA, nrow=length(x), ncol=length(levels(x)) )
	dimnames(result)[[2]] <- paste(".",levels(x),sep="")
	for (i in 1:length(levels(x))) {
		result[,i] <- as.integer( x==levels(x)[i] )
	}
	return(result)
}

bracket <- function( x0, y0, x1=x0, y1=y0, 
						offset=1, length=offset/2, 
						side=1, col="grey", ... ) {
	switch( side, 
		{ # side=1 below
			segments( x0=x0, y0=y0+length-offset, x1=x0, y1=y0-offset, col=col, ...)
			segments( x0=x0, y0=y0-offset, x1=x1, y1=y0-offset, col=col, ...)
			segments( x0=x1, y0=y1-offset, x1=x1, y1=y0+length-offset, col=col, ...)
			},
		{ # side=2 left
			segments( x0=x0+length-offset, y0=y0, x1=x0-offset, y1=y0, col=col, ...)
			segments( x0=x0-offset, y0=y0, x1=x1-offset, y1=y1, col=col, ...)
			segments( x0=x0-offset, y0=y1, x1=x0+length-offset, y1=y1, col=col, ...)
			},
		{ # side=3 above
			segments( x0=x0, y0=y0-length+offset, x1=x0, y1=y0+offset, col=col, ...)
			segments( x0=x0, y0=y0+offset, x1=x1, y1=y0+offset, col=col, ...)
			segments( x0=x1, y0=y1+offset, x1=x1, y1=y0-length+offset, col=col, ...)
			},
		{ # side=4 right
			segments( x0=x0-length+offset, y0=y0, x1=x0+offset, y1=y0, col=col, ...)
			segments( x0=x0+offset, y0=y0, x1=x1+offset, y1=y1, col=col, ...)
			segments( x0=x0+offset, y0=y1, x1=x0-length+offset, y1=y1, col=col, ...)
			} ) # end of switch
}

is.inrange <- function ( x, range=c(0,1) ) {    
   return( range[1]<=x & x<=range[2] ) 
}

### CONVERSION FUNCTIONS ###

f2st <- function( hz, base=50 ) {
   semi1 <- log( 2^(1/12) ) # equals log(2)/12
   return( ( log(hz) - log(base) ) / semi1 )
}

st2f <- function( st, base=50 ) {
   semi1 <- log( 2^(1/12) ) # equals log(2)/12
   return( exp(st*semi1) * base )
}

f2bark <- function(hz) { 26.81 / (1+(1960/hz) ) -0.53 }
bark2f <- function(bark)  { 1960* (bark+0.53) / (26.28-bark) }

f2mel <- function (hz) { 2595 * log(1 + hz/700, base = 10) }
mel2f <- function (mel) { 700 * (10^(mel/2595) - 1) }
