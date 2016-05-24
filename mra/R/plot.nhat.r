plot.nhat <- function( x, ci=TRUE, smooth=TRUE, occasions=-1, smubass=5,...){

	# Plot N hat
	y <- x$n.hat
	occasion <- cumsum( c(1,x$intervals ) )   # seq( along=y ) this plots all intervals as equal length
	ns <- length( y )
	if( length(names(y)) != ns ) names(y) <- 1:ns
	nms <- names(y)
	if( any(occasions <= 0) ){
		occasions <- 1:length(y)
	}
	occasions <- occasions[ occasions != 1 ]

	y <- y[occasions]
	occasion <- occasion[occasions]
	nms   <- nms[occasions]

	n.hat <- y
	if( ci ){
		lower.ci <- x$n.hat.lower[occasions]
		upper.ci <- x$n.hat.upper[occasions]

		plot( range(occasion), range(n.hat,lower.ci,upper.ci),
			type="n", xaxt="n", xlab="Occasion", ylab="N estimate", ...)
		axis( side=1, at=occasion, labels=nms )
		lines(occasion, lower.ci, type="l", lty=2)
		lines(occasion, upper.ci, type="l", lty=2)
		lines(occasion, n.hat, type="b")
	} else {

		plot(occasion, n.hat, type="b", xaxt="n",  ...)
		axis(1, at=occasion, labels=nms )

	}


	if( smooth & exists("supsmu") ){
		sm <- supsmu( occasion, n.hat, bass=smubass )
		lines( sm, lwd=3, col=3 )
	} else {
		sm <- NA
	}

	invisible(sm)
}
