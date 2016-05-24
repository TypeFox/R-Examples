## typically arise as output from "benchmark" function
## objects of class "turbosim" consist of a named list with the following components:
	## method.names = vector of unique identifiers for the methods being compared
	## fail = matrix of logical values where (i,j)-entry indicates whether algorithm j of simulation iteration i failed (produced an error)
	## conv = matrix of logical values where (i,j)-entry indicates whether algorithm j of simulation iteration i converged (i.e. did not exceed maximum number of iterations or runtime)
	## value.objfn = matrix of minimizers of the objective function
	## runtime = matrix of computation times
	## method = vector of algoritm names = c("em", "squarem", "pem", "decme", "qn")
	## control.method = list of control parameters used in the turboem function of the simulation
	## control.run = list of run parameters used in the turboem function of the simulation

########################################
########################################
## print ###############################
########################################
########################################

print.turbosim <- function(x, ...) {
	cat("\nBenchmark study over", nrow(x$fail), "repetitions.\n")
	cat("\nMethods:\n")
	for(i in seq_along(x$method.names)) {
		cat(i, ". ", x$method.names[i], "\n", sep="")
	}
	cat("\nFunctions to summarize and visualize results:\n    summary(), boxplot(), dataprof(), pairs()\n")
}

########################################
########################################
## table of failures ###################
########################################
########################################

summary.turbosim <- function(object, which.methods = seq_along(object$method), method.names = object$method.names[which.methods], eps = 0.1, sol = NULL, ...) {
	x <- object
	
	## sol = true objective function value, if known
	nmethods <- length(method.names)
	value.objfn <- x$value.objfn[,which.methods]
	fail <- x$fail[,which.methods]
	conv <- x$conv[,which.methods]
	if(length(method.names) != length(which.methods)) {
		stop("length of method.names must be equal to length of which.methods")
	}
	
	if(nmethods > 1) {
	
		if(!is.null(sol)) {
			fail3 <- value.objfn > eps + sol
		} else if(is.null(sol)) {
			fail3 <- value.objfn > eps + replicate(nmethods, abs(apply(value.objfn, 1, min, na.rm=TRUE)))
		}
		tab <- rbind(colSums(fail), colSums(!conv, na.rm=TRUE), colSums(fail3, na.rm=TRUE))
		
		if(x$control.run$stoptype=="maxiter") {
			name2 <- paste("Exceeded", x$control.run$maxiter, "iter. ")
		} else if(x$control.run$stoptype=="maxtime") {
			name2 <- paste("Exceeded", round(x$control.run$maxtime/60,2), "min. ")
		} else if(x$control.run$stoptype=="user") {
			name2 <- paste("Did not converge")
		} 
		
		rownames(tab) <- c("Algorithm failed ", name2, paste("objfn > min(objfn) +", eps, ""))
		colnames(tab) <- method.names
		return(t(tab))
	} else {
		df.method <- data.frame(fail=fail, convergence=conv, value.objfn=value.objfn, itr=x$itr[,which.methods], fpeval=x$fpeval[,which.methods], objfeval=x$objfeval[,which.methods], had.error=!x$errors[,which.methods] %in% c("",NA))
		return(summary(df.method))
	}
}	

########################################
########################################
## boxplots of runtimes ################
########################################
########################################

boxplot.turbosim <- function(x, which.methods = seq_along(x$method), method.names = x$method.names[which.methods], whichfail = (x$fail | !x$conv)[,which.methods], xunit="sec", log=FALSE, ...) {
	if(ncol(whichfail) > length(which.methods)) {
		whichfail <- whichfail[,which.methods]
	}
	if(length(method.names) != length(which.methods)) {
		stop("length of method.names must be equal to length of which.methods")
	}

	runtime <- x$runtime[,which.methods]
	## Modified (JFB Feb2012)
	fail <- whichfail
	runtime <- ifelse(fail, NA, runtime)
	# fail <- x$fail[,which.methods]
	opar <- par(mar=c(4.1, 13.1, 2.1, 2.1))
	on.exit(par(opar))
	if(xunit=="sec") {
		xlab <- "time (sec)"
		dat <- runtime
	} else if(xunit=="min") {
		xlab <- "time (min)"
		dat <- runtime/60
	}
	if(log) {
		dat <- log(dat)
		xlab <- paste("log", xlab)
	}
	boxplot(dat, horizontal=TRUE, axes=FALSE, main="Runtimes ", xlab=xlab)
	axis(1)
	axis(2, las=1, at=1:length(method.names), labels=paste(method.names, " (", colSums(!fail), ")", sep=""))
	mtext("Algorithm (# iter)", side=2, line=11.5, cex=1.5)
}

########################################
########################################
## data profile ########################
########################################
########################################

dataprof <- function(x, ...) {
	UseMethod("dataprof")
}
dataprof.turbosim <- function(x, which.methods = seq_along(x$method), method.names = x$method.names[which.methods], whichfail = (x$fail | !x$conv)[,which.methods], col, lty, nout = 50, xlim, ...) {
	## whichfail = matrix of the same dimension of 'runtimes' identifying which were failures
	## nout = number of points over which the empirical distribution function is calculated
	##require(RColorBrewer)
	method <- x$method[which.methods]
	nmethod <- length(which.methods)
	fail <- x$fail[,which.methods]
	runtime <- x$runtime[,which.methods]
	if(ncol(whichfail) > length(which.methods)) {
		whichfail <- whichfail[,which.methods]
	}
	if(length(method.names) != length(which.methods)) {
		stop("length of method.names must be equal to length of which.methods")
	}

	nsim <- nrow(fail)
	unq <- unique(method)
	if(missing(col)) {
		##col <- 1:nmethod
		col <- unlist(sapply(1:length(unq), function(u) rep(u, table(method)[unq][u])))
	}
	if(missing(lty)) {
		##lty <- 1:nmethod
		lty <- unlist(sapply(1:length(unq), function(u) 1:table(method)[unq][u]))
	}
	xs <- freq <- matrix(NA, nout, nmethod)
	for(k in 1:nmethod) {
		xs[,k] <- quantile(runtime[,k], probs=seq(0, 1, length.out=nout), na.rm=TRUE)
		freq[,k] <- sapply(xs[,k], function(u) sum(!whichfail[,k] & runtime[,k] <= u))
	}
	xs <- rbind(xs, rep(max(runtime, na.rm=TRUE), nmethod))
	freq <- rbind(freq, freq[nrow(freq),])
	if(missing(xlim)) {
		xlim <- range(runtime, na.rm=TRUE)/60
	}
	titl <- "Distribution of T = 'Time to Convergence' \ngiven no failures"
	##opar <- par(mar=c(5.1, 4.1, 4.1, 12.1))
	opar <- par(mar=c(5.1, 4.1, 4.1, 2.1 + max(nchar(method.names))-1))
	on.exit(par(opar))
	plot(freq[,1], type="n", xlim=xlim, ylim=range(freq/nsim), xlab="Time t (minutes)", ylab="Proportion of iterations with T < t", main=titl, ...)
	for(k in 1:nmethod) {
		lines(xs[,k]/60, freq[,k]/nsim, col=col[k], lty=lty[k], lwd=2)
	}
	if(x$control.run$stoptype=="maxtime") segments(4/5*x$control.run$maxtime/60, 1, x$control.run$maxtime/60, 1, col="gray")
	legend(par("usr")[2], par("usr")[4], paste(method.names, " (", colSums(!whichfail, na.rm=TRUE), ")", sep=""), col=col, lty=lty, bty="n", lwd=2, xpd=TRUE)
	##legend("bottomright", paste(method.names, " (", colSums(!whichfail, na.rm=TRUE), ")", sep=""), col=col, lty=lty, bty="n", lwd=2)
}

########################################
########################################
## scatterplot matrix ##################
########################################
########################################

pairs.turbosim <- function(x, which.methods=seq_along(x$method), method.names = x$method.names[which.methods], whichfail = (x$fail | !x$conv)[,which.methods], ...) {
	## which.methods = vector identifying which columns of 'runtimes' will be used, e.g. which of the which.methods will be included in the scatterplot matrix
	## whichfail = matrix of the same dimension of 'runtimes' identifying which were failures

	if(length(method.names) != length(which.methods)) {
		stop("length of method.names must be equal to length of which.methods")
	}

	runtimes <- x$runtime[,which.methods]
	nmeth <- ncol(runtimes)
	if(ncol(whichfail)==length(x$method.names)) {
		whichfail <- whichfail[,which.methods]
	}

	opar <- par(mfcol=c(nmeth,nmeth+1), mar=c(1,0.5,0.5,0.5))
	on.exit(par(opar))
	legend.loc <- ceiling(length(which.methods)/2)
	for(k in 1:nmeth) {
		plot(1, 1, type="n", axes=FALSE, xlab="", ylab="")
		if(k == ceiling(length(which.methods)/2)) {
			legend("center", c("y < x", "y > x", "L1 fit"), col=c("red", "blue", "black"), pch=c(19,19,NA), lwd=c(NA,NA,2), cex=1.5, bty="n")
		}
	}
	for(j in 1:nmeth) {
		for(i in 1:nmeth) {
			if(i != j) {
				par(mar=c(4.1, 4.1, 2.1, 1.1))
				failj <- whichfail[,j]==1
				faili <- whichfail[,i]==1
				x <- runtimes[,j]
				y <- runtimes[,i]
				sel <- !faili | !failj
				if(any(sel)) {
					x[failj] <- ifelse(all(failj), max(y[!faili], na.rm=TRUE)+1, max(x[!failj], na.rm=TRUE)) ##max(x[!failj], na.rm=TRUE)
					y[faili] <- ifelse(all(faili), max(x[!failj], na.rm=TRUE)+1, max(y[!faili], na.rm=TRUE))
					type <- "p"
					##reg <- try(lm(y ~ 0 + x))
					reg <- try(rq(y ~ 0 + x, tau=0.5))
					##reg <- try(rq(y ~ x, tau=0.5))
					cols <- ifelse(y > x, "blue", "red")
					cols[failj] <- "red"
					cols[faili] <- "blue"
					x <- x[sel]
					y <- y[sel]
					cols <- cols[sel]
					plot(y ~ x, xlab="", ylab="", main="", col=cols, pch=19, type=type, axes=FALSE, ...)
					axis(1)
					axis(2)
					title(xlab=method.names[j], ylab=method.names[i], line=2.5)
					##abline(0,1, col="purple", lwd=2)
					if(class(reg) != "try-error") {
						title(main=paste("y = ", round(coef(reg)["x"], 2), "x", sep=""), line=1)
						##title(main=paste("y = ", round(coef(reg)["(Intercept)"], 2), " + ", round(coef(reg)["x"], 2), "x", sep=""), line=1)
						abline(reg=reg, col="black", lwd=2)
					}
					##lines(lowess(x, y, f=0.5), col="brown", lwd=2)
				} else {
					plot(1, 1, type="n", axes=FALSE, xlab="", ylab="")
					text(1, 1, "all failures")
				}
			} else if(i==j) {
				par(mar=c(1,0.5,0.5,0.5))
				plot(1, 1, type="n", axes=FALSE, xlab="", ylab="")
				text(1, 1, method.names[i], cex=18/(nchar(method.names[i])+6))
				text(1, 0.85, paste("(", sum(!whichfail[,i]), ")", sep=""), cex=18/(nchar(method.names[i])+7))
			}
		}
	}
	# par(mar=c(5.1, 4.1, 4.1, 2.1))
	# par(mfcol=c(1,1))
}
