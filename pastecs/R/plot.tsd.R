"plot.tsd" <-
function (x, series=1, stack=TRUE, resid=TRUE, col=par("col"), lty=par("lty"), labels=dimnames(X)[[2]], leg = TRUE, lpos=c(0,0), xlab="time", ylab="series", main=paste("Series decomposition by", x$specs$method, "-", x$specs$type), ...) {
    ser <- x$ts
    if (length(ser) == 1) {		# We have the decomposition of a single series in the object
    	if (series != 1)
    		stop("the series does not exist in this object")
    	sers <- unclass(x$series)
    } else {					# We have the decomposition of several series in the object
    	if (series < 1 || series > length(ser))
    		stop("the series does not exist in this object")
    	sers <- unclass(x$series[[series]])
    }
    model.type <- x$specs$type		# Additive or multiplicative
    if (is.null(model.type))
    	model.type <- "additive"	# By default
    ncomp <- ncol(sers)
    if (model.type == "additive") {
    	series <- drop(sers %*% rep(1, ncomp))
    	# Equivalent to: series <- apply(sers, 1, sum)
    } else {
    	series <- apply(sers, 1, prod)
    }
    X <- cbind(series = series, sers)
    # Does the last series represent residuals?
    if (dimnames(sers)[[2]][ncomp] == "residuals") {
    	if (resid == TRUE) nplot <- ncomp + 1 else nplot <- ncomp
    } else {		# No residuals
    	nplot <- ncomp + 1
    	resid <- FALSE
    }
    col <- rep(col, nplot)
    lty <- rep(lty, nplot)
    if (stack == TRUE) {
    	oldpar <- par("mar", "oma", "mfrow", "tck")
    	on.exit(par(oldpar))
    	par(mar = c(0, 6, 0, 6), oma = c(6, 0, 4, 0), tck = -0.01)
    	par(mfrow = c(nplot, 1))
    	for (i in 1:nplot) {
    		plot.type <- if (i < nplot | resid == FALSE) "l" else if (model.type == "additive") "h" else "p"
    		
    		# Solution using segments for residuals based on 1
    		#x.1 <- rnorm(1000,10)  # random numbers around 10
			#plot(x.1,type='n')
			#segments(seq(x.1),10,seq(x.1),x.1)  # plot segments basing on 10
    		
        	plot(X[, i], type = plot.type, xlab = "", ylab = "", axes = FALSE, col = col[i], lty = lty[i], ...)
        	if (i == nplot & resid == TRUE & min(X[, i]) < 0 & max(X[, i]) > 0) 
            	if (model.type == "additive") abline(h = 0) else abline(h = 1)
            box()
        	right <- i%%2 == 0
        	axis(2, labels = !right)
        	axis(4, labels = right)
        	mtext(labels[i], 2, 3)
        }
    	axis(1, labels = TRUE)
    	axis(3, labels = FALSE)
    	mtext(xlab, 1, 3)
    } else {		# Stack is false
		X <- as.ts(X)
    	if (resid == TRUE) {
			X <- as.ts(X)
		} else {
			n <- ncol(X) - 1	
			X <- as.ts(X[,1:n])	
		}
		ts.plot(X, gpars=list(col=col, lty=lty, xlab=xlab, ylab=ylab, main=main))
    	if (leg == TRUE) legend(lpos[1], lpos[2], labels[1:n], col=col[1:n], lty=lty[1:n])
    }
    invisible()
}
