"lines.escouf" <-
function(x, level=x$level, lhorz=TRUE, lvert=TRUE, lvars=TRUE, col=2, lty=2, ...) {
	# The next function actually draw the lines
	escouf.lines <- function(X, Level, Lhorz, Lvert, Lvars, Col, Lty, ...) {
		if (is.null(Level)) { 	# Missing level argument
			stop("You must provide a value for level!")
		} else {				# We draw the lines
			n <- length(X$RV)
			if (Lhorz==TRUE)		# We draw also an horizontal line
				lines(c(1, n), c(Level, Level), lty=Lty, col=Col)
			# How many variables do we keep?
			nvars <- length(X$RV[X$RV<Level])
			if (Lvert==TRUE)		# We draw also a vertical line
				lines(c(nvars+0.5, nvars+0.5), c(-0.1,1.5), lty=Lty, col=Col)
			if (Lvars==TRUE)		# We change also colors of selected variables labels
				axis(1, 1:nvars, labels=as.character(X$vr[1:nvars]), col.axis=Col)
		}
	}
	invisible(escouf.lines(x, level, lhorz, lvert, lvars, col, lty, ...))
}
