getTimestamp <- function() {
	.Call("getTimestamp", package="datautils")
}

getElapsed <- function(stamp) {
	.Call("getElapsed", stamp, package="datautils")
}


# build a "merged" list,
# ie 2-level list tree :
# - 1st level has as many cells as each parameter,
# - 2nd level is filled with a list of the i-th elements of each parameter. 
mergeToList <- function(...) {
	thecall <- as.list(sys.call())
	len <- length(thecall)
	res <- list()
	refsize <- NA
	if (len == 1) return (res)
	
	# function that fills the list to return
	appendfunc <- function(arg) {
		arg <- eval(arg)
		if (is.na(refsize)) {
			refsize <<- length(arg)
			for (I in 1:length(arg)) res[[I]] <<- list() # <<- modifies enclosing env.
		}	
		if (length(arg) != refsize) stop("arguments' lengths must match")

		for (I in 1:length(arg)) {
			if (is.list(arg)) {
				elt <- arg[[I]]
			} else {
				elt <- arg[I]
			}
			res[[I]] <<- c(res[[I]], list(elt))
		} 
	}
	
	# process all arguments
	lapply(thecall[2:len], appendfunc)
	return(res)
}


# return the purity of ground truth labels wrt inferred clusters
getPurity <- function(truthLabels, inferLabels) {
	if (length(truthLabels) != length(inferLabels)) stop("truthLabels and inferLabels lengths")
	getMajorityCount <- function(vec) {
		# given a (factor-castable) vector of labels, return the count of cells that equals the most frequently present value.
		vec <- as.factor(vec)
		lev <- levels(vec)
		maxcount <- -Inf
		for (mind in 1:length(lev)) {
			if (length(which(vec==lev[mind])) > maxcount) {
				maxcount <- length(which(vec==lev[mind]))
			}
		}
		return(maxcount)
	}	
		
		
	inferlevels <- levels(as.factor(inferLabels))
	labgroups <- lapply(inferlevels, function(x) which(inferLabels==x))
	reslist <- lapply(labgroups, function(x) getMajorityCount(truthLabels[x]))
	res <- numeric()
	for(mind in 1:length(reslist)) res[mind] <- reslist[[mind]]
	return(sum(res) / length(truthLabels))
}

# hack of the plotmeans function of gplots package, to allow native scale if the grouping variable is numeric.
# introduction of xvals and nummeans parameter.
plotmeanshack <- function (formula, data = NULL, subset, na.action, bars = TRUE, 
    p = 0.95, minsd = 0, minbar = NULL, maxbar = NULL, xlab = names(mf)[2], 
    ylab = names(mf)[1], mean.labels = FALSE, ci.label = FALSE, 
    n.label = TRUE, digits = getOption("digits"), col = "black", 
    barwidth = 1, barcol = "blue", connect = TRUE, ccol = col, 
    legends = names(means), xaxt, use.t = TRUE, nummeans=TRUE, ...) 
{

    is.R <- get("is.R")
    if (is.null(is.R)) 
        is.R <- function(x) FALSE
    if (!is.R()) {
        if (col == "black") 
            col <- 1
        if (barcol == "blue") 
            barcol <- 2
    }
    if (invalid(formula) || (length(formula) != 3)) 
        stop("formula missing or incorrect")
    if (invalid(na.action)) 
        na.action <- options("na.action")
    m <- match.call(expand.dots = FALSE)
    if (is.R()) {
        if (is.matrix(eval(m$data, parent.frame()))) 
            m$data <- as.data.frame(data)
    }
    else {
        if (is.matrix(eval(m$data, FALSE))) 
            m$data <- as.data.frame(data)
    }
    m$... <- m$bars <- m$barcol <- m$p <- NULL
    m$minsd <- m$minbar <- m$maxbar <- NULL
    m$xlab <- m$ylab <- NULL
    m$col <- m$barwidth <- NULL
    m$digits <- m$mean.labels <- m$ci.label <- m$n.label <- NULL
    m$connect <- m$ccol <- m$legends <- m$labels <- NULL
    m$xaxt <- m$use.t <- m$nummeans <- NULL
    m[[1]] <- as.name("model.frame")
    mf <- eval(m, parent.frame())
    response <- attr(attr(mf, "terms"), "response")
    wFact <- which(attr(attr(mf, "terms"), "dataClasses") == 
        "factor")
    for (i in wFact) mf[, i] <- factor(mf[, i])
    means <- sapply(split(mf[[response]], mf[[-response]]), mean, 
        na.rm = TRUE)

    if (nummeans) {
    	xvals <- as.numeric(names(means))
    } else {
    	xvals <- 1:length(means)
    } 
    ns <- sapply(sapply(split(mf[[response]], mf[[-response]]), 
        na.omit, simplify = FALSE), length)

    xlim <- c(xvals[1] - 0.5, xvals[length(xvals)] + 0.5)
    if (!bars) {
        plot(x=xvals, y=means, ..., col = col, xlim = xlim)
    }
    else {
        myvar <- function(x) var(x[!is.na(x)])
        vars <- sapply(split(mf[[response]], mf[[-response]]), 
            myvar)
        vars <- ifelse(vars < (minsd^2), (minsd^2), vars)
        if (use.t) 
            ci.width <- qt((1 + p)/2, ns - 1) * sqrt(vars/ns)
        else ci.width <- qnorm((1 + p)/2) * sqrt(vars/ns)
        if (length(mean.labels) == 1) {
            if (mean.labels == TRUE) 
                mean.labels <- format(round(means, digits = digits))
            else if (mean.labels == FALSE) 
                mean.labels <- NULL
        }

        
        plotCI(x = xvals, y = means, uiw = ci.width, 
            xaxt = "n", xlab = xlab, ylab = ylab, labels = mean.labels, 
            col = col, xlim = xlim, lwd = barwidth, barcol = barcol, 
            minbar = minbar, maxbar = maxbar, ...)
        if (invalid(xaxt) || xaxt != "n") 
            axis(1, at = xvals, labels = legends)
        if (ci.label) {
            ci.lower <- means - ci.width
            ci.upper <- means + ci.width
            if (!invalid(minbar)) 
                ci.lower <- ifelse(ci.lower < minbar, minbar, 
                  ci.lower)
            if (!invalid(maxbar)) 
                ci.upper <- ifelse(ci.upper > maxbar, maxbar, 
                  ci.upper)
            labels.lower <- paste(" \n", format(round(ci.lower, 
                digits = digits)), sep = "")
            labels.upper <- paste(format(round(ci.upper, digits = digits)), 
                "\n ", sep = "")
            text(x = xvals, y = ci.lower, labels = labels.lower, 
                col = col)
            text(x = xvals, y = ci.upper, labels = labels.upper, 
                col = col)
        }
    }
    if (n.label) 
        if (is.R()) 
            text(x = xvals, y = par("usr")[3], labels = paste("n=", 
                ns, "\n", sep = ""))
        else {
            axisadj <- (par("usr")[4] - (par("usr")[3]))/75
            text(x = xvals, y = par("usr")[3] + axisadj, 
                labels = paste("n=", ns, "\n", sep = ""))
        }
    if (!invalid(connect) & !identical(connect, FALSE)) {
        if (is.list(connect)) {
            if (length(ccol) == 1) 
                ccol <- rep(ccol, length(connect))
            for (which in 1:length(connect)) lines(x = connect[[which]], 
                y = means[connect[[which]]], col = ccol[which])
        }
        else lines(x=xvals, y=means, ..., col = ccol)
    }
}




plot.deldir <- function (x, add = FALSE, wlines = c("both", "triang", "tess"), 
    wpoints = c("both", "real", "dummy", "none"), number = FALSE, 
    cex = 1, nex = 1, col = NULL, lty = NULL, pch = NULL, xlim = NULL, 
    ylim = NULL, xlab = "x", ylab = "y", showrect = FALSE, fill = NULL, ...) 
{
    if (!inherits(x, "deldir")) 
        stop("Argument \"x\" is not of class deldir.\n")
    wlines <- match.arg(wlines)
    wpoints <- match.arg(wpoints)
    col <- if (is.null(col)) 
        c(1, 1, 1, 1, 1)
    else rep(col, length.out = 5)
    lty <- if (is.null(lty)) 
        1:2
    else rep(lty, length.out = 2)
    pch <- if (is.null(pch)) 
        1:2
    else rep(pch, length.out = 2)
    plot.del <- switch(wlines, both = TRUE, triang = TRUE, tess = FALSE) # good example of how to use switch
    plot.dir <- switch(wlines, both = TRUE, triang = FALSE, tess = TRUE)
    plot.rl <- switch(wpoints, both = TRUE, real = TRUE, dummy = FALSE, 
        none = FALSE)
    plot.dum <- switch(wpoints, both = TRUE, real = FALSE, dummy = TRUE, 
        none = FALSE)
    delsgs <- x$delsgs
    dirsgs <- x$dirsgs
    n <- x$n.data
    rw <- x$rw
    if (plot.del) {
        x1 <- delsgs[, 1]
        y1 <- delsgs[, 2]
        x2 <- delsgs[, 3]
        y2 <- delsgs[, 4]
    }
    if (plot.dir) {
        u1 <- dirsgs[, 1]
        v1 <- dirsgs[, 2]
        u2 <- dirsgs[, 3]
        v2 <- dirsgs[, 4]
    }
    X <- x$summary[, "x"]
    Y <- x$summary[, "y"]
    if (!add) {
        pty.save <- par()$pty
        on.exit(par(pty = pty.save))
        par(pty = "s")
        if (is.null(xlim)) 
            xlim <- rw[1:2]
        if (is.null(ylim)) 
            ylim <- rw[3:4]
        plot(0, 0, type = "n", xlim = xlim, ylim = ylim, xlab = xlab, 
            ylab = ylab, axes = FALSE, ...)
        axis(side = 1)
        axis(side = 2)
    }
    if (plot.del) 
        segments(x1, y1, x2, y2, col = col[1], lty = lty[1], 
            ...)
    if (plot.dir) {
		# handle all polygons and colors in one single polygon() command.
		# if fill is NULL, set with white
		
		umin <- min(c(u1, u2))
		vmin <- min(c(v1, v2))
		umax <- max(c(u1, u2))
		vmax <- max(c(v1, v2))
		
		if (is.null(fill)) fill <- rgb(1,1,1) # use rgb format
		polygons <- matrix(nrow=0, ncol=2)
		for(i in 1:n) {
			inds <- which(dirsgs$ind1 == i) # get all tesselation segments matching a single point
			inds <- unique(c(inds, which(dirsgs$ind2 == i)))
			tmppts <- rbind(cbind(u1[inds], v1[inds]), cbind(u2[inds], v2[inds]))
			# manage the "corner" case : if a min and a max are attained by points in tmppts on both dimension, 
			# add the "corner" to the convex hull.
			if (any(tmppts[,1] == umin) && any(tmppts[,2] == vmin)) tmppts <- rbind(tmppts, c(umin, vmin))
			if (any(tmppts[,1] == umin) && any(tmppts[,2] == vmax)) tmppts <- rbind(tmppts, c(umin, vmax))
			if (any(tmppts[,1] == umax) && any(tmppts[,2] == vmin)) tmppts <- rbind(tmppts, c(umax, vmin))
			if (any(tmppts[,1] == umax) && any(tmppts[,2] == vmax)) tmppts <- rbind(tmppts, c(umax, vmax))
			 
        	inds <- chull(tmppts) # convex hull in clockwise order
        	tmppts <- tmppts[inds,]	
			polygons <- rbind(polygons, tmppts, rep(NA,2)) # close each polygon with NA values
		} # warning : for does not define a lexical scope. 
		polygon(polygons[,1], polygons[,2], col=fill)
	}
    if (plot.rl) {
        x.real <- X[1:n]
        y.real <- Y[1:n]
        points(x.real, y.real, pch = pch[1], col = col[3], cex = cex, 
            ...)
    }
    if (plot.dum) {
        x.dumm <- X[-(1:n)]
        y.dumm <- Y[-(1:n)]
        points(x.dumm, y.dumm, pch = pch[2], col = col[4], cex = cex, 
            ...)
    }
    if (number) {
        xoff <- 0.02 * diff(range(X))
        yoff <- 0.02 * diff(range(Y))
        text(X + xoff, Y + yoff, 1:length(X), cex = nex, col = col[5], 
            ...)
    }
    if (showrect) 
        do.call(rect, as.list(x$rw)[c(1, 3, 2, 4)])
    invisible()
}


upper <- function(d) {
	if (d < 2) {
		stop("d should at least equal 2")
	}
	res <- matrix(nrow=0, ncol=2)
	for (i in 1:(d-1)) {
		cur <- seq(i+1, d)
		res <- rbind(res, cbind(rep(i, length(cur)), cur))
	}
	return(res)
}	

	
	
	





