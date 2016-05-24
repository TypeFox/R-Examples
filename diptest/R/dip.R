### S-interface to Hartigan's algorithm for "The dip test for unimodality"
###
### Beginning:	Dario Ringach <dario@wotan.cns.nyu.edu>
### Rest:	Martin Maechler <maechler@stat.math.ethz.ch>

dip <- function(x, full.result = FALSE, min.is.0 = FALSE, debug = FALSE)
{
    allRes <- (!is.logical(rFull <- full.result))
    if(allRes) {
	if(full.result %in% c("all"))
	    rFull <- TRUE
	else stop(gettextf("'full.result' = \"%s\"' is not valid", full.result))
    }
    if(rFull) cl <- match.call()

    if(is.unsorted(x))
	x <- sort(x, method="quick")
    n <- as.integer(length(x))
    r <- .C(diptst,
	    x	= as.double(x),
	    n	= n,
	    dip = double(1),
	    lo.hi = integer(4),
	    ifault= integer(1),
	    gcm =   integer(n),
	    lcm =   integer(n),
	    mn	=   integer(n),
	    mj	=   integer(n),
	    min.is.0 = as.logical(min.is.0),
	    debug = as.integer(debug)# FALSE/TRUE or 2, 3, ...
            )[if(rFull) TRUE else "dip"]
    if(rFull) {
	l.GL <- r$lo.hi[3:4]
	length(r$gcm) <- l.GL[1]
	length(r$lcm) <- l.GL[2]
	length(r$lo.hi) <- 2L
	u <- x[r$lo.hi]
	structure(class = "dip",
		  c(list(call = cl), r,
		    list(xl = u[1], xu = u[2], full.result=full.result),
		    if(allRes) getCM(r$mn, r$mj, n)))
    }
    else r[[1]]
}

getCM <- function(mn, mj, n = length(mn)) {
    stopifnot(length(mn) <= n, length(mj) <= n) # currently '=='...
    ## First recover "the full GCM / LCM" - by repeating what happened in C
    ## in the first "loop" :
    low <- 1L ; high <- n
    gcm <- lcm <- integer(n) # pre-allocate!  {maybe smaller ?}

    ## Collect the change points for the GCM from HIGH to LOW. */
    gcm[i <- 1L] <- high
    while(gcm[i] > low)
	gcm[(i <- i+1L)] <- mn[gcm[i]]
    length(gcm) <- i

    ## Collect the change points for the LCM from LOW to HIGH. */
    lcm[i <- 1L] <- low
    while(lcm[i] < high)
	lcm[(i <- i+1L)] <- mj[lcm[i]]
    length(lcm) <- i
    list(GCM = gcm, LCM = lcm)
}

print.dip <- function(x, digits = getOption("digits"), ...)
{
    stopifnot(is.integer(n <- x$n), is.numeric(D <- x$dip),
	      length(lh <- x$lo.hi) == 2)
    cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
	"\n\n", sep = "")
    xLU.c <- sapply(x$x[lh], formatC, digits=digits, width=1)
    cat("n = ", n,".  Dip statistic, D_n = ",
	format(D, digits=digits)," = ",
	format(2*n* D, digits=digits),"/(2n)\n",
	sprintf(" Modal interval [xL, xU] = [x[%d], x[%d]] = [%s, %s]\n",
		lh[1], lh[2], xLU.c[1], xLU.c[2]),
	sprintf(" GCM and LCM have %d and %d nodes inside [xL, xU], respectively",
                ## 3 5 7 9 1 3 5 7
		length(x$gcm), length(x$lcm)),
	if(x$full.result == "all")
	sprintf(", and\n%17s %d and %d nodes in   [x_1, x_n].\n", "",
		length(x$GCM), length(x$LCM)) else ".\n",
	sep="")
    invisible(x)
}

aLine <- function(r.dip, lType = c("gcm","lcm","GCM","LCM"),
		  type = "b", col="red3", lwd=1.5, ...)
{
    lType <- match.arg(lType)
    stopifnot(is.numeric(x <- r.dip$x), is.integer(r.dip$n),
	      is.integer(i <- r.dip[[lType]]) # 'gcm' or 'lcm' or component
	      )
    e <- if(lType %in% c("gcm","GCM")) .01*min(diff(unique(x))) else 0
    i <- i[i != 0]
    lines(x[i], ecdf(x)(x[i] - e),
	  type=type, col=col, lwd=lwd, ...)
}

plot.dip <- function(x, do.points=(n < 20), ## <- plot.stepfun()
		     colG="red3", colL="blue3", colM="forest green",
		     col.points=par("col"), col.hor=col.points, ## <- plot.stepfun():
		     doModal=TRUE, doLegend=TRUE, ...)
{
    stopifnot(is.integer(n <- x$n), is.numeric(D <- x$dip),
	      length(lh <- x$lo.hi) == 2)
    Fn <- ecdf(x$x)
    ## and now manipulate the call such that it's plotted nicely
    cl <- x$call[1:2]
    cl[[1]] <- as.name("ecdf") ; names(cl)[2] <- ""
    attr(Fn, "call") <- cl
    chD <- formatC(D, digits=pmax(3, getOption("digits")-2))
    tit <- bquote("Dip" ~~ {D[n] == D[.(n)]} == .(chD))
    plot(Fn, do.points=do.points, col.points=col.points, col.hor=col.hor,
	 verticals=TRUE, col.vert = "sky blue", lwd=2, ...)
    title(tit, adj = 0, line = 1.25)
    aLine(x, "gcm", col=colG)
    aLine(x, "lcm", col=colL)
    if(doCM.2 <- (x$full.result == "all")) {
        aLine(x, "GCM", col=colG, lty=5)
        aLine(x, "LCM", col=colL, lty=5)
    }
    if(doModal) {
	x12 <- x$x[lh]
	abline(v= x12, col = colM, lty = 2)
	op <- par(mgp = c(3, 1/16, 0))# should not need on.exit(par(op)) here ..
	axis(3, at=x12, labels = expression(x[L], x[U]),
	     tick=FALSE, col.axis = colM)
	par(op)
    }
    if(doLegend) {
	txt <- c("greatest convex minorant GCM",
### make sure have *no* [TAB] in next string !
		 "least     concave majorant LCM")
	t1 <- paste(txt," in [xL, xU]")
	legend("topleft", bty = "n",
	       if(doCM.2) c(t1, txt) else t1,
	       lwd=1.5, col = c(colG, colL), lty= if(doCM.2) c(1,1,5,5) else 1)
    }
    invisible()
}
