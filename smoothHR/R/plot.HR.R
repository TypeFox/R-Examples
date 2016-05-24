plot.HR <- function(x, predictor, prob=NULL, pred.value=NULL, conf.level=0.95, round.x=NULL, ref.label=NULL, col, main, xlab, ylab, lty, xlim, ylim, xx, ...) {
	object <- x
	if ( !inherits(object, "HR") ) stop("Object must be of class HR")
	mydata <- object$dataset
	fit <- object$coxfit
	if ( missing(round.x) ) round.x <- 5
	if ( !missing(pred.value) ) prob <- 0.5
	if ( !missing(prob) ) if(prob < 0 | prob > 1) stop("The argument 'prob' must be between o and 1")
	if ( missing(prob) & missing(pred.value) ) prob <- 0
	if ( !missing(pred.value) & !missing(xlim) ) if ( pred.value < min(xlim) | pred.value > max(xlim) ) stop("The reference value is out of range of 'xlim'")
	if ( missing(predictor) ) stop("Missing predictor")
	if ( missing(col) ) col <- c("black", "black", "grey85")
	if ( missing(ylab) ) ylab <- c("Ln HR(Z,Zref)")
	if ( missing(lty) ) lty <- c(1,3)
	ctype <- "FALSE"
	qvalue <- (1+conf.level)/2
	linear.predictor <- FALSE
	k1 <- 9999
	k <- which(names(mydata) == predictor)
	k <- c(k, k1)
	if (k[1] == 9999) stop ("predictor must be in data")
	k <- k[1]
	a <- mydata
	#if (!missing(pred.value) & (pred.value<min(a[,k]) | pred.value>max(a[,k]))) stop("The reference value is out of range of x)
	if ( missing(xlab) ) xlab <- names(a)[k]
	n.predictor <- names(a)[k]
	n <- dim(a)[1]
	if (prob == 0) {
		eta.no.ref <- predict(fit,type = "terms")
		if ( inherits(eta.no.ref, "numeric") ) {kp <- 1; eta.no.ref <- cbind(eta.no.ref,eta.no.ref);}
		else {kp <- grep( predictor, colnames(eta.no.ref) );}
		eta.xref <- min(eta.no.ref[,kp])
		ii <- which.min(eta.no.ref[,kp])
		xref <- a[ii,k]
		eta.ref <- eta.no.ref[,kp]-eta.xref
		indices <- grep(names(a)[k], dimnames(fit$x)[[2]])
		submatriz.diseno <- fit$x[,indices]
		if (is.matrix(submatriz.diseno) == FALSE) linear.predictor <- TRUE
		submatriz.var <- fit$var[indices, indices]
		xref1 <- rep(fit$x[ii,indices], dim(fit$x)[1])
		if (linear.predictor == FALSE) {
			xref1 <- matrix(xref1, nrow=dim(fit$x)[1], ncol=dim(submatriz.diseno)[2], byrow=TRUE)
		}
		if (linear.predictor == TRUE) {
			xref1 <- matrix(xref1, nrow=dim(fit$x)[1], ncol=1, byrow=TRUE)
		}
		eta.ref1 <- fit$x[,indices]-xref1
		var.eta.ref1 <- rep(NA, n)
		for (i in 1:n) var.eta.ref1[i] <- eta.ref1[i,]%*%fit$var[indices,indices]%*%eta.ref1[i,]
		se.eta.ref1 <- sqrt(var.eta.ref1)
	}
	if (prob > 0 & prob < 1) {
		eta.no.ref <- predict(fit, type="terms")
		if ( inherits(eta.no.ref, "numeric") ) {kp <- 1; eta.no.ref <- cbind(eta.no.ref, eta.no.ref);}
		else {kp <- grep( predictor, colnames(eta.no.ref) )}
		ord <- order(a[,k])
		if ( !missing(pred.value) ) {
			pp <- seq(0, 1, len=1000)
			app <- quantile(a[,k], pp)
			qq <- which(app<=pred.value)
			qq1 <- max(qq)
			prob <- qq1/1000
		}
		ind.prob <- trunc(prob*n)
		xref <- a[,k][ord[ind.prob]]
		eta.xref <- eta.no.ref[,kp][ord[ind.prob]]
		eta.ref <- eta.no.ref[,kp]-eta.xref
		indices <- grep(names(a)[k], dimnames(fit$x)[[2]])
		submatriz.diseno <- fit$x[,indices]
		if (is.matrix(submatriz.diseno) == FALSE) linear.predictor <- TRUE
		submatriz.var <- fit$var[indices, indices]
		xref1 <- rep(fit$x[ord[ind.prob],indices], dim(fit$x)[1])
		if (linear.predictor == FALSE) {
			xref1 <- matrix(xref1, nrow=dim(fit$x)[1], ncol=dim(submatriz.diseno)[2], byrow=TRUE)
		}
		if (linear.predictor == TRUE) {
			xref1 <- matrix(xref1, nrow=dim(fit$x)[1], ncol=1, byrow=TRUE)
		}
		eta.ref1 <- fit$x[,indices]-xref1
		var.eta.ref1 <- rep(NA,n)
		for (i in 1:n) var.eta.ref1[i] <- eta.ref1[i,]%*%fit$var[indices,indices]%*%eta.ref1[i,]
		se.eta.ref1 <- sqrt(var.eta.ref1)
	}
	if (prob == 1) {
		eta.no.ref <- predict(fit, type="terms")
		if ( inherits(eta.no.ref, "numeric") ) {kp <- 1; eta.no.ref <- cbind(eta.no.ref, eta.no.ref);}
		else {kp <- grep( predictor, colnames(eta.no.ref) )}
		eta.xref <- max(eta.no.ref[,kp])
		ii <- which.max(eta.no.ref[,kp])
		xref <- a[ii,k]
		eta.ref <- eta.no.ref[,kp]-eta.xref
		indices <- grep(names(a)[k], dimnames(fit$x)[[2]])
		submatriz.diseno <- fit$x[,indices]
		if (is.matrix(submatriz.diseno) == FALSE) linear.predictor <- TRUE
		submatriz.var <- fit$var[indices,indices]
		xref1 <- rep(fit$x[ii,indices], dim(fit$x)[1])
		if (linear.predictor == FALSE) {
			xref1 <- matrix(xref1, nrow=dim(fit$x)[1], ncol=dim(submatriz.diseno)[2], byrow=TRUE)
		}
		if (linear.predictor == TRUE) {
			xref1 <- matrix(xref1, nrow=dim(fit$x)[1], ncol=1, byrow=TRUE)
		}
		eta.ref1 <- fit$x[,indices]-xref1
		var.eta.ref1 <- rep(NA,n)
		for (i in 1:n) var.eta.ref1[i] <- eta.ref1[i,]%*%fit$var[indices,indices]%*%eta.ref1[i,]
		se.eta.ref1 <- sqrt(var.eta.ref1)
	}
	if ( missing(main) ) main <- paste("Smooth log hazard ratio for", names(a)[k])
	tmat <- cbind(eta.ref, eta.ref-qnorm(qvalue)*se.eta.ref1, eta.ref+qnorm(qvalue)*se.eta.ref1)
	line <- rep(0, n)
	jj <- match(sort(unique(a[,k])), a[,k])
	if ( missing(xlim) ) xlim <- c( min(a[,k]), max(a[,k]) )
	else {
		if ( missing(ylim) ) {
			index1 <- which( a[jj,k] >= min(xlim) & a[jj,k] <= max(xlim) )
			index <- jj[index1]
			ylim <- c( min(tmat[index,2]), max(tmat[index,3]) )
		}
	}
	if ( missing(ylim) ) ylim <- c( min(tmat[,2]), max(tmat[,3]) )
	if ( xref < min(a[,k]) | xref > max(a[,k]) ) stop("The reference value is out of range of x")
	if ( xref < min(xlim) | xref > max(xlim) ) stop("The reference value is out of range of 'xlim'")
	matplot(a[jj,k], tmat[jj,], type="l", lty=c(1, 5, 5, 2), xaxt="n", ylim=ylim, xlim=xlim, xlab=xlab, ylab=ylab, col=c(1, 2, 2, 1), main=main, ...)
	xxx <- round( seq(min(a[,k]), max(a[,k]),len=5) )
	if ( missing(xx) ) xx <- c( min(a[,k]), round(xref,1), xxx[2], xxx[3], xxx[4], max(a[,k]) )
	axis(1, xx, ...)
	m <- length(jj)
	x <- rep(NA, 2*m+1)
	y <- rep(NA, 2*m+1)
	for (l in 1:m) {
		x[l] <- a[jj,k][l]
		x[m+l] <- a[jj,k][m+1-l]
		y[l] <- tmat[jj,2][l]
		y[m+l] <- tmat[jj,1][m+1-l]
	}
	x[m+1] <- x[m]
	x[2*m+1] <- x[1]
	y[2*m+1] <- tmat[jj,2][1]
	polygon(c(x), c(y), col=col[3], ...)
	y <- rep(NA, 2*m+1)
	for (l in 1:m) {
		x[l] <- a[jj,k][l]
		x[m+l] <- a[jj,k][m+1-l]
		y[l] <- tmat[jj,3][l]
		y[m+l] <- tmat[jj,1][m+1-l]
	}
	x[m+1] <- x[m]
	x[2*m+1] <- x[1]
	y[2*m+1] <- tmat[jj,2][1]
	polygon(c(x), c(y), col=col[3], ...)
	#_____________________________
	y <- rep(NA, 2*m+1)
	for (l in 1:m) {
		x[l] <- a[jj,k][l]
		x[m+l] <- a[jj,k][m+1-l]
		y[l] <- tmat[jj,3][l]
		y[m+l] <- tmat[jj,2][m+1-l]
	}
	x[m+1] <- x[m]
	x[2*m+1] <- x[1]
	y[2*m+1] <- tmat[jj,2][1]
	polygon(c(x), c(y), col=col[3], border="white", ...)
	#_____________________________
	y <- rep(NA, 2*m+1)
	for (l in 1:m) {
		x[l] <- a[jj,k][l]
		x[m+l] <- a[jj,k][m+1-l]
		y[l] <- tmat[jj,3][l]
		y[m+l] <- tmat[jj,2][m+1-l]
	}
	x[m+1] <- x[m]
	x[2*m+1] <- x[1]
	y[2*m+1] <- tmat[jj,2][1]
	polygon(c(x), c(y), col=col[3], border=col[2], lty=lty[2], lwd=1.5, ...)
	# (lty=3 poiwise confidence bands; lty=4 : -.-)
	#_____________________________
	x <- rep(NA, 2*m+1)
	y <- rep(NA, 2*m+1)
	for (l in 1:m) {
		x[l] <- a[jj,k][l]
		x[m+l] <- a[jj,k][m+1-l]
		y[l] <- tmat[jj,1][l]
		y[m+l] <- tmat[jj,1][m+1-l]
	}
	x[m+1] <- x[m]
	x[2*m+1] <- x[1]
	y[2*m+1] <- tmat[jj,1][1]
	polygon(c(x), c(y), col=col[3], border=col[1], lty=lty[1], ...)
	#_______________________________
	abline(0, 0, lty=2)
	abline(v=min(a[,k]), col="white")
	abline(v=max(a[,k]), col="white")
	if ( missing(xlim) ) {
		v1 <- min(a[,k])+( max(a[,k])-min(a[,k]) )/10
		v2 <- min(a[,k])+9*( max(a[,k])-min(a[,k]) )/10
	} else {
		v1 <- min(xlim)+( max(xlim)-min(xlim) )/10
		v2 <- min(xlim)+9*( max(xlim)-min(xlim) )/10
	}
	if ( missing(ylim) ) {
		y[1] <- max(tmat[,3])/2
		y[2] <- min(tmat[,2])
	} else {
		y[1] <- max(ylim)/2
		y[2] <- min(ylim)
	}
	if ( !missing(ref.label) ) n.predictor <- ref.label
	if (xref > v1 & xref < v2) {
		arrows(xref, y[1], xref, y[2], length=0.08)
		ys <- y[1]
		if (ys > 2*y[1]-(2*y[1]-y[2])/10)  {
			text(xref, y[1], paste(n.predictor, "=", round( xref, round.x) ), adj=c(0.5, 2.3), ...)
		}
		if (ys <= 2*y[1]-(2*y[1]- y[2])/10)  {
			text(xref, y[1], paste( n.predictor, "=", round(xref,round.x) ), adj=c(0.5, -0.7), ...)
		}
	}
	if (xref <= v1) {
		v3 <- ( max(xlim)-min(xlim) )/100
		xref2 <- xref
		if ( xref == min(xlim) ) xref2 <- xref+min(0.05, v3)
		arrows(xref2, y[1], xref2, y[2], length=0.08)
		ys <- y[1]
		if (ys > 2*y[1]-(2*y[1]-y[2])/10)  {
			text(xref, y[1], paste( n.predictor, "=", round(xref,round.x) ), adj=c(0, 2.3), ...)
		}
		if (ys <= 2*y[1]-(2*y[1]-y[2])/10) {
			text(xref, y[1], paste( n.predictor, "=", round(xref,round.x) ), adj=c(0, -0.7), ...)
		}
	}
	if (xref >= v2) {
		v3 <- ( max(xlim)-min(xlim) )/100
		xref2 <- xref
		if ( xref == max(xlim) ) xref2 <- xref-min(0.05, v3)
		arrows(xref2, y[1], xref2, y[2], length=0.08)
		ys <- y[1]
		if (ys > 2*y[1]-(2*y[1]-y[2])/10) {
			text(xref, y[1], paste( n.predictor, "=", round(xref,round.x) ), adj=c(1, 2.3), ...)
		}
		if (ys <= 2*y[1]-(2*y[1]-y[2])/10) {
			text(xref, y[1], paste( n.predictor, "=", round(xref,round.x) ), adj=c(1, -0.7), ...)
		}
	}
}
