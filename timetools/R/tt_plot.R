points.TimeIntervalDataFrame <- function (x, y=NULL, cursor=NULL, type='p',
					  lty=1:6, lwd=1, pch=1:25, col=NULL, ...) {
	if (is.null(cursor)) {
		if (is.null(y)) y <- names(x)
		y <- data.frame (x[y])

		# individual serial parameters
		lty <- rep (lty, length.out=length(y))
		lwd <- rep (lwd, length.out=length(y))
		if( is.data.frame(pch) )
			pch <- as.list(pch[rep(1:length(pch),
					       length.out=length(y))]) else
			pch <- rep (pch, length.out=length(y))
		if (is.null(col))
			col <- sample (colors(), length(y)) else
			col <- rep (col, length.out=length(y))
		
		# plotting
		if (type %in% c('l', 'b', 'o') ) {
			trash <- mapply (segments, y0=y[-nrow(x),,drop=FALSE],
				y1=y[-1,,drop=FALSE], col=col, lty=lty, lwd=lwd,
				MoreArgs=list(x1=start(x)[-1], x0=end(x)[-nrow(x)]) )
		}
		trash <- mapply (segments, y0=y, col=col, lty=lty, lwd=lwd,
			MoreArgs=list(x0=start(x), x1=end(x)) )
		if (type %in% c('p', 'b', 'o', 'h')) {
			type <- ifelse (type %in% c('b', 'o'), 'p', type)
			trash <- mapply (points, y=y, pch=pch, col=col, ...,
				MoreArgs=list(x=start(x), type=type))
			trash <- mapply (points, y=y, pch=pch, col=col, ...,
				MoreArgs=list(x=end(x), type=type))
		}
		pch <- sapply(pch, paste, collapse=', ')
		if (length (y) > 0)
			invisible(data.frame (names=names(y), type, lty, lwd, pch,
					col, stringsAsFactors=FALSE))
	} else {
		x <- as.TimeInstantDataFrame (x, cursor=cursor)
		invisible(points (x=x, y=y, type=type, lty=lty, lwd=lwd, pch=pch,
				col=col, ...))
	}
}

lines.TimeIntervalDataFrame <- function (x, y=NULL, cursor=NULL, type='l',
					 lty=1:6, lwd=1, pch=1:25, col=NULL, ...)
{
	invisible (points (x=x, y=y, cursor=cursor, type=type, lty=lty, lwd=lwd,
			   pch=pch, col=col, ...))
}

plot.TimeIntervalDataFrame <- function (x, y=NULL, cursor=NULL,
					type='p', lty=1:6, lwd=1, pch=1:25, col=NULL,
					xlim=NULL, ylim=NULL,
					log='', main='', sub='', xlab='', ylab='',
					ann=par('ann'), axes=TRUE, asp=NA, ...) {
	if (is.null(cursor))
	{
		if (is.null(y)) y <- names(x)
		y <- data.frame (x[y])

		# overal graph parameters
		if (is.null(xlim)) xlim <- range(c(start(x), end(x)))
		if (is.null(ylim)) ylim <- compute.lim (unlist(y), na.rm=TRUE)

		# prepare plot
		plot (NA, xlim=xlim, ylim=ylim, log=log, main=main, sub=sub,
		      xlab=xlab, ylab=ylab, ann=ann, axes=FALSE, asp=asp, ...)
		if (axes) {
			axis(2)
			axis.POSIXct (1, c(start(x), end(x)) ) 
			box()
		}
		invisible (points (x=x, y=names(y), cursor=cursor, type=type,
				   lty=lty, lwd=lwd, pch=pch, col=col, ...))
	} else {
		x <- as.TimeInstantDataFrame (x, cursor=cursor)
		invisible(plot (x=x, y=y, type=type, lty=lty, lwd=lwd, pch=pch,
			col=col,
  			xlim=xlim, ylim=ylim, log=log,
  			main=main, sub=sub, xlab=xlab, ylab=ylab, ann=ann,
  			axes=axes, asp=asp, ...))
	}
}

barplot.TimeIntervalDataFrame <- function(height, format='', ...) {
	dots <- list(...)
	if( !'names.arg' %in% names(dots) )
		dots$names.arg <- paste(
			format(start(height), format=format),
			format(end(height), format=format),
			sep='\n')
	height <- t(as.data.frame(height))
	do.call(barplot, c(height=list(height), dots))
}

points.TimeInstantDataFrame <- function (x, y=NULL, type='p', lty=1:6,
					 lwd=1, pch=1:25, col=NULL, ...) {
	if (is.null(y)) y <- names(x)
	y <- data.frame (x[y])

	# individual serial parameters
	type <- rep (type, length.out=length(y))
	lty <- rep (lty, length.out=length(y))
	lwd <- rep (lwd, length.out=length(y))
	if( is.data.frame(pch) )
		pch <- as.list(pch[rep(1:length(pch), length.out=length(y))]) else
		pch <- rep (pch, length.out=length(y))
	if (is.null(col))
		col <- sample (colors(), length(y)) else
		col <- rep (col, length.out=length(y))
	
	# plotting
	trash <- mapply (points, y=y, pch=pch, col=col, type=type, lty=lty,
			 lwd=lwd, ...,
			MoreArgs=list(x=when(x)))
	pch <- sapply(pch, paste, collapse=', ')
	if (length(y) > 0)
		invisible(data.frame (names=names(y), type, lty, lwd, pch, col,
				stringsAsFactors=FALSE))
}

lines.TimeInstantDataFrame <- function (x, y=NULL, type='l', lty=1:6, lwd=1,
					pch=1:25, col=NULL, ...) {
	invisible (points (x=x, y=y, type=type, lty=lty, lwd=lwd, pch=pch,
			   col=col, ...))
}

plot.TimeInstantDataFrame <- function (x, y=NULL,
					type='p', lty=1:6, lwd=1, pch=1:25,
					col=NULL,
					xlim=NULL, ylim=NULL,
					log='', main='', sub='', xlab='', ylab='',
					ann=par('ann'), axes=TRUE, asp=NA, ...)
{
	if (is.null(y)) y <- names(x)
	y <- data.frame (x[y])

	# overal graph parameters
	if (is.null(xlim)) xlim <- range(when(x))
	if (is.null(ylim)) ylim <- compute.lim (unlist(y), na.rm=TRUE)

	# prepare plot
	plot (NA, xlim=xlim, ylim=ylim, log=log, main=main, sub=sub,
	      xlab=xlab, ylab=ylab, ann=ann, axes=FALSE, asp=asp, ...)
	if (axes) {
		axis(2)
		axis.POSIXct (1, when(x)) 
		box()
	}
	invisible (points (x=x, y=names(y), type=type, lty=lty, lwd=lwd, pch=pch,
			   col=col, ...))
}

barplot.TimeInstantDataFrame <- function(height, format='', ...) {
	dots <- list(...)
	if( !'names.arg' %in% names(dots) )
		dots$names.arg <- format(when(height), format=format)
	height <- t(as.data.frame(height))
	do.call(barplot, c(height=list(height), dots))
}

points.SubtimeDataFrame <- function (x, y=NULL,
					type='p', lty=1:6, lwd=1, pch=1:25,
					col=NULL,
					as.is=TRUE, ...) {
	if (is.null(y)) y <- names(x)
	y <- data.frame (x[y])

	if( as.is ) abscisses <- 1:nrow(x) else
		abscisses <- as.numeric(when(x))

	# individual serial parameters
	type <- rep (type, length.out=length(y))
	lty <- rep (lty, length.out=length(y))
	lwd <- rep (lwd, length.out=length(y))
	if( is.data.frame(pch) )
		pch <- as.list(pch[rep(1:length(pch), length.out=length(y))]) else
		pch <- rep (pch, length.out=length(y))
	if (is.null(col))
		col <- sample (colors(), length(y)) else
		col <- rep (col, length.out=length(y))
	
	# plotting
	trash <- mapply (points, y=y, pch=pch, col=col, type=type, lwd=lwd,
			 lty=lty, ...,
 			 MoreArgs=list(x=abscisses))
	pch <- sapply(pch, paste, collapse=', ')
	if (length(y) > 0)
		invisible(data.frame (names=names(y), type, lty, lwd, pch, col,
				stringsAsFactors=FALSE))
}

lines.SubtimeDataFrame <- function (x, y=NULL,
					type='l', lty=1:6, lwd=1, pch=1:25,
					col=NULL,
					as.is=TRUE, ...) {
	invisible (points (x=x, y=y, type=type, lty=lty, lwd=lwd, pch=pch, col=col,
			   as.is=as.is, ...))
}

plot.SubtimeDataFrame <- function (x, y=NULL,
					type='p', lty=1:6, lwd=1, pch=1:25,
					col=NULL,
					xlim=NULL, ylim=NULL,
					log='', main='', sub='', xlab='', ylab='',
					ann=par('ann'), axes=TRUE, asp=NA,
					as.is=TRUE, format=NULL, ...)
{
	if (is.null(y)) y <- names(x)
	y <- data.frame (x[y])

	if( as.is ) abscisses <- 1:nrow(x) else
		abscisses <- as.numeric(when(x))

	# overal graph parameters
	if (is.null(xlim)) xlim <- range(abscisses)
	if (is.null(ylim)) ylim <- compute.lim (unlist(y), na.rm=TRUE)

	# prepare plot
	plot (NA, xlim=xlim, ylim=ylim, log=log, main=main, sub=sub,
	      xlab=xlab, ylab=ylab, ann=ann, axes=FALSE, asp=asp, ...)
	if (axes) {
		axis(2)
		axis(1, abscisses, labels=format(when(x), format))
		box()
	}
	invisible (points (x=x, y=names(y), type=type, lty=lty, lwd=lwd, pch=pch,
			   col=col, as.is=as.is, ...))
}

barplot.SubtimeDataFrame <- function(height, format=NULL, ...) {
	dots <- list(...)
	if( !'names.arg' %in% names(dots) )
		dots$names.arg <- format(when(height), format=format)
	height <- t(as.data.frame(height))
	do.call(barplot, c(height=list(height), dots))
}

