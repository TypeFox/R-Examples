"plot.escouf" <-
function(x, level=x$level, lhorz=TRUE, lvert=TRUE, lvars=TRUE, lcol=2, llty=2, diff=TRUE, dlab="RV' (units not shown)", dcol=4, dlty=par("lty"), dpos=0.8, type="s", xlab="variables", ylab="RV", main=paste("Escoufier's equivalent vectors for:",x$data), ...) {
	# The next function actually draw the graph
	escouf.graph <- function(X, Level, Lhorz, Lvert, Lvars, Lcol, Llty, Diff, Dlab, Dcol, Dlty, Dpos, Type, Xlab, Ylab, Main, ...) {
		n <- length(X$RV)
		plot(X$RV, type=Type, xlab=Xlab, ylab=Ylab, main=Main, xaxs="i", xaxt="n", ...)
		axis(1, 1:length(X$RV), labels=as.character(X$vr))
		# Do we plot also RV.diff?
		if (Diff==TRUE) {
			# Calculate RV.diff
			RV.diff <- X$RV[2:n] - X$RV[1:(n-1)]
			# Scale RV.diff to the same range as RV
			RVd <- RV.diff
			RVds <- (RVd-min(RVd))/max(RVd)*(max(X$RV)-min(X$RV))+min(X$RV)
			# Plot the line
			lines((1:(n-1))+0.5, RVds, col=Dcol, lty=Dlty)
			# Draw the label
			xPos <- n*Dpos
			yInd <- round(xPos); if (yInd<length(RVds)) yInd <- length(RVds)
			text(xPos, RVds[yInd], Dlab, pos=3, col=Dcol)
		}
		if (is.null(Level)==FALSE) { # We draw the lines
			if (Lhorz==TRUE)		# We draw also a horizontal line
				lines(c(1, n), c(Level, Level), lty=Llty, col=Lcol)
			# How many variables do we keep?
			nvars <- length(X$RV[X$RV<Level])
			if (Lvert==TRUE)		# We draw also a vertical line
				lines(c(nvars+0.5, nvars+0.5), c(-0.1,1.5), lty=Llty, col=Lcol)
			if (Lvars==TRUE)		# We change also colors of selected variables labels
				axis(1, 1:nvars, labels=as.character(X$vr[1:nvars]), col.axis=Lcol)
		}
			
	}
	invisible(escouf.graph(x, level, lhorz, lvert, lvars, lcol, llty, diff, dlab, dcol, dlty, dpos, type, xlab, ylab, main, ...))
}
