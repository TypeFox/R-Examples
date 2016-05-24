# Copyright (C) 2011 Pierrick Bruneau, see README for full notice


mymvn2plot <- function (w, mu, sigma, k = 15, alone = FALSE, col = NA, alphacol=0.8, alphanocol=0.5, lty="solid") {
    p <- length(mu)
    if (p != 2) 
        stop("two-dimensional case only")
    if (any(unique(dim(sigma)) != p)) 
        stop("mu and sigma are incompatible")
    ev <- eigen(sigma, symmetric = TRUE)
    s <- sqrt(rev(sort(ev$values)))
    V <- t(ev$vectors[, rev(order(ev$values))])
    theta <- (0:k) * (pi/(2 * k))
    x <- s[1] * cos(theta)
    y <- s[2] * sin(theta)
    xy <- cbind(c(x, -x, -x, x), c(y, y, -y, -y))
	# OK : ici les vecteurs V sont en ligne, donc on prend cos(t) * v1 + sin(t) * v2
	# ce qui correspond bien a la courbe parametrique d'un cercle/ellipse.
    xy <- xy %*% V
    xy <- sweep(xy, MARGIN = 2, STATS = mu, FUN = "+")
    if (alone) {
        xymin <- apply(xy, 2, FUN = "min")
        xymax <- apply(xy, 2, FUN = "max")
        r <- ceiling(max(xymax - xymin)/2)
        xymid <- (xymin + xymax)/2
        plot(xy[, 1], xy[, 2], xlim = c(-r, r) + xymid[1], ylim = c(-r, 
            r) + xymid[2], xlab = "x", ylab = "y", type = "n")
    }

    l <- length(x)
    i <- (l+1):(2*l)

	xy[i,] <- apply(xy[i,], 2, rev)
	
	i <- (3*l+1):(4*l)
	xy[i,] <- apply(xy[i,], 2, rev)
		
	i <- 1:(4*l)
	# rescale w in order to keep grey components
	w <- 1 - w
	
	# if a specific color is provided (text format), mix grey weight specified by w with this color.
	if(!is.na(col)) {
		w <- w * 255
		col <- col2rgb(col)
		col <- (col + w) / 2
		col <- col / 255
		alpha <- alphacol
	} else {
		col <- matrix(rep(w, 3))
		alpha <- alphanocol
	}
	polygon(xy[,1], xy[,2], border="black", lty=lty, col=rgb(col[1,1],col[2,1],col[3,1], alpha))

	# set xy order so as to define a path

    x <- s[1]
    y <- s[2]
    xy <- cbind(c(x, -x, 0, 0), c(0, 0, y, -y))
    xy <- xy %*% V
    xy <- sweep(xy, MARGIN = 2, STATS = mu, FUN = "+")
    invisible()
}



mySmoothScatter <- function(data, model, xlim, ylim) {
	smoothScatter(data, nrpoints=0, xlim=xlim, ylim=ylim)
	if(!is.null(model)) {
		for(i in 1:length(model$w)) {
			mymvn2plot(model$w[i], model$mean[[i]], model$cov[[i]])
		}
	}
}

getColor <- function(index) {
	# color index set selected from colors() returned values
	# 1st color=black
	cols=c(24,74,32,26, 57,31, 83, 101, 90, 116, 114, 126, 551)
	cycle <- sum((index %/% 13) > 0) # hack modulo for a regular cycle: 1..i..n 1..i...
	
	return(colors()[cols[(index %% 13) + cycle]])
}	

displayScatter <- function(data=NULL, model=NULL, labels=NULL, datasizes=NULL, compcolors=NULL, complabels=NULL, compstrokes="solid", space=1:2, xlim=NULL, ylim=NULL, main="", xlab="", ylab="", smooth=FALSE, alphacol=0.8, alphanocol=0.5, cex.lab=1, lwd=1) {
	# parameters:
	# - data: a data set in nxd matrix to be plot. Cast to 2D (two first)
	# - model: a GMM to be plot. Cast to 2D (two first)
	# - labels: an optional numeric label set for data elements. May be reset if inconsistent size. indexes are related to colors returned by getColor()
	# - datasizes: if scalar=>set cex size of all elements. if vector=> pch for each element. May be reset if inconsistent size.
	# - compcolors: if scalar=> color to associate with component regions. May be reset if inconsistent size. indexes are related to colors returned by getColor()
	# - complabels: text labels to associate with components. May be reset.
	# - xlim, ylim: range for plot. if not specified, will be infered from data.
	# - main, xlab & ylab: transmitted as such to plot function
	
	par(oma=c(0,0,0,0), mar=c(4.5, 5.5, 1.5, 1.5), mex=0.5)
	
	# either model or data should be not null.
	if(is.null(data) && is.null(model)) {
		stop("either data or model should be not null")
	}	
	
	# check model first. Set to void if NULL.
	if(is.null(model)) model <- newGmm()
	
	# model should be a GMM
	if(is.null(model$w) || is.null(model$mean) || is.null(model$cov)) {
		stop("model should be a GMM")
	}
		
	# sizes must be consistent
	# set k.
	k1 <- length(model$w)
	k2 <- length(model$mean)
	k3 <- length(model$cov)
	if((k1 != k2) || (k1 != k3)) {
		stop("inconsistency in model size")
	}
	k <- k1

	# cast means before going further
	if(k>0) {
		for(i in 1:k) {
			model$mean[[i]] <- as.numeric(model$mean[[i]])
		}
	}


	
	# set default data set if void
	if(is.null(data)) {
		n <- k
		data <- matrix(model$mean[[1]], nrow=1, ncol=length(model$mean[[1]]))
		if(k>1) {
			for(i in 2:k) {
				data <- rbind(data, model$mean[[i]])
			}
		}
	} else {
		n <- dim(data)[1]
	}


	# check space argument
	if(length(space) != 2) stop("length of space argument should be 2")

	# cast data and model to 2D, according to space argument
	data <- data[,space]
	if(k>0) {
		for(i in 1:k) {
			model$mean[[i]] <- model$mean[[i]][space]
			model$cov[[i]] <- model$cov[[i]][space, space]
		}
	}

	# detect xlim and ylim if null
	if(is.null(xlim)) {
		xlim[1] <- min(data[,1])
		xlim[2] <- max(data[,1])
		#margin <- (xlim[2] - xlim[1]) * 0.1
		#xlim[1] <- xlim[1] - margin
		#xlim[2] <- xlim[2] + margin
	}
	
	if(is.null(ylim)) {
		ylim[1] <- min(data[,2])
		ylim[2] <- max(data[,2])
		#margin <- (ylim[2] - ylim[1]) * 0.1
		#ylim[1] <- ylim[1] - margin
		#ylim[2] <- ylim[2] + margin

	}
	
	# get labels from model or specific param.
	if(is.null(labels)) {
		labels <- model$labels
	}
	
	# fill if null, check size
	if(is.null(labels)) {
		# all points are black
		labels <- rep(1, n)
	} else {
		# avoid black colors for points
		labels <- labels + 1
	}
	if(length(labels) != n) stop('inconsistent length for labels')
	
	# replace numeric labels (=color indexes) with names
	newlabels <- character()
	for(i in 1:n) {
		newlabels[i] <- getColor(labels[i])
	}
	labels <- newlabels
	
	# datasizes may be null. If not: either a scalar or vector. If a vector, should be of size n.
	if(is.null(datasizes)) datasizes <- 1
	if(length(datasizes) == 1)  datasizes <- rep(datasizes, n)
	if(length(datasizes) != n) stop("inappropriate datasizes argument")
	
	# compcolors may be null. If null, set with a void numeric().
	if(is.null(compcolors)) compcolors <- numeric()
	if(length(compcolors) == 1) compcolors <- rep(compcolors, k)
	if((length(compcolors) != k) && (length(compcolors) != 0)) stop("inappropriate compcolors argument")
	
	if(length(compstrokes) == 1) compstrokes <- rep(compstrokes, k)
	if(length(compstrokes) != k) stop("inappropriate compstrokes argument")

	
	# replace component colors with names
	# manage NA values
	newcompcolors <- character()
	if(length(compcolors) > 0) {
		newcompcolors[k] <- NA
		for(i in 1:k) {
			if(!is.na(compcolors[i])) newcompcolors[i] <- getColor(compcolors[i])
		}
		compcolors <- newcompcolors
	}
	
	# check complabels
	if(is.null(complabels))  complabels <- ""
	if(length(complabels) == 1) complabels <- rep(complabels, k)
	if(length(complabels) != k) stop("inappropriate complabels argument")
	
	# 2 types of plots: classic or smoothed. labels and cex are not used for smooth.
	if(!smooth) {
		plot(data, xlim=xlim, ylim=ylim, pch=4, col=labels, cex=datasizes, main=main, xlab=xlab, ylab=ylab, cex.lab=cex.lab, lwd=lwd)
	} else {
		smoothScatter(data, nrpoints=0, xlim=xlim, ylim=ylim, main=main, xlab=xlab, ylab=ylab, cex.lab=cex.lab, bandwidth=0.3)
	}
	
	# hack : 2 successive loops for printing opaque (colored) components first
	if(k>0) {
		for(i in 1:k) {
			if(!is.na(compcolors[i])) {
				mymvn2plot(model$w[i], mu=model$mean[[i]], sigma=model$cov[[i]], col=compcolors[i], alphacol=alphacol, lty=compstrokes[i])
				text(complabels[i], x=model$mean[[i]][1], y=model$mean[[i]][2])
			}
		}
		for(i in 1:k) {
			if(is.na(compcolors[i])) {
				mymvn2plot(model$w[i], mu=model$mean[[i]], sigma=model$cov[[i]], col=compcolors[i], alphanocol=alphanocol, lty=compstrokes[i])
				text(complabels[i], x=model$mean[[i]][1], y=model$mean[[i]][2])
			}
		}
	}
}






gridGen <- function(xlim=c(-10, 10), ylim=c(-10, 10), step=50) {
	# generates 2D matrix of coordinates of regular grid of elements
	stepsizeX <- (xlim[2] - xlim[1]) / step
	stepsizeY <- (ylim[2] - ylim[1]) / step
	res <- matrix(nrow=step^2, ncol=2)
	for(i in 0:(step-1)) {
		for(j in 0:(step-1)) {
			# manage flipping coordinates between screen and matrix convention
			res[i*step + j + 1,] <- c(xlim[1] + (i+1/2)*stepsizeX, ylim[2] - (j+1/2)*stepsizeY)
		}
	}
	return(res)
}

# mettre des factors pour variables categorielle, sinon SVM fait une espece de regression...
# build data frame from usual format
buildFrame <- function(datamatrix, labels, dims=1:2) {
	datamatrix <- as.data.frame(datamatrix)
	labels <- as.factor(as.numeric(labels))
	datamatrix <- datamatrix[,dims]
	datamatrix <- cbind(datamatrix, labels)
	return(datamatrix)
}

displaySVM <- function(svm.model, dataframe, displayPoints=TRUE, subset=NULL, steps=100, alpha=0.4, lwd=1) {
	# restricted to 2D cases
	# set limits
	xlim <- ylim <- numeric()
	xlim[1] <- min(dataframe[,1])
	xlim[2] <- max(dataframe[,1])
	ylim[1] <- min(dataframe[,2])
	ylim[2] <- max(dataframe[,2])
	
	xstep <- (xlim[2] - xlim[1]) / steps
	ystep <- (ylim[2] - ylim[1]) / steps
	

	if(is.null(subset)) subset <- 1:dim(dataframe)[1]

	names <- names(dataframe)

	par(oma=c(0,0,0,0), mar=c(4.5, 5.5, 1.5, 1.5), mex=0.5, cex.lab=1.5)	
	plot(type="n", xlim=xlim, ylim=ylim, x=NULL, y=NULL, xlab=names[1], ylab=names[2])

	# build grid data to predict
	datagrid <- gridGen(xlim, ylim, steps)
	gridlabs <- rep(1, dim(datagrid)[1])
	datagrid <- buildFrame(datagrid, gridlabs)
	names(datagrid) <- names
	
	svm.pred <- predict(svm.model, datagrid)
	svm.pred <- as.numeric(svm.pred)
	# avoid black
	svm.pred <- svm.pred + 1

	# set up grid viewport
	pos <- par("plt")
	grid::pushViewport(grid::viewport(x=grid::unit(pos[1], "npc"), y=grid::unit(pos[3], "npc"), 
		width=grid::unit(pos[2]-pos[1], "npc"), height=grid::unit(pos[4]-pos[3], "npc"), 
		just=c("left", "bottom")))
	grid::pushViewport(grid::viewport(xscale=xlim, yscale=ylim))


	# create alpha blended colors from labels
	# OPTIM: simplify loop
	#gridcolors <- numeric()
	#for(i in 1:dim(datagrid)[1]) {	
	#	curcolor <- getColor(svm.pred[i])
	#	curcolor <- as.numeric(col2rgb(curcolor)) / 255
	#	curcolor <- rgb(curcolor[1], curcolor[2], curcolor[3], alpha)
	#	gridcolors[i] <- curcolor
	#}
	tempcolors <- getColor(svm.pred)
	tempcolors <- col2rgb(tempcolors) / 255
	gridcolors <- rgb(tempcolors[1,], tempcolors[2,], tempcolors[3,], alpha)
	
	
	# plot decision regions
	grid::grid.rect(x=grid::unit(datagrid[,1], "native"), y=grid::unit(datagrid[,2], "native"), 
		width=grid::unit(xstep, "native"), height=grid::unit(ystep, "native"), 
		gp=grid::gpar(fill=gridcolors, lty=0, col=gridcolors))
		
	# then plot points, if not null
	if(displayPoints) {
		labs <- as.numeric(dataframe[subset,3]) + 1
		newlabs <- character()
		for(i in 1:length(labs)) {
			newlabs[i] <- getColor(labs[i])
		}
		labs <- newlabs
		points(dataframe[subset,1], dataframe[subset,2], col=labs, pch=4, lwd=lwd)
	}
}


displayNnet <- function(nnet.model, datamatrix, datalabels, subset=NULL, displayPoints=TRUE, steps=100, alpha=0.4, lwd=1) {	
	# restricted to 2D cases
	# set limits
	
	xlim <- ylim <- numeric()
	xlim[1] <- min(datamatrix[,1])
	xlim[2] <- max(datamatrix[,1])
	ylim[1] <- min(datamatrix[,2])
	ylim[2] <- max(datamatrix[,2])
	
	xstep <- (xlim[2] - xlim[1]) / steps
	ystep <- (ylim[2] - ylim[1]) / steps

	if(is.null(subset)) subset <- 1:dim(datamatrix)[1]

	par(oma=c(0,0,0,0), mar=c(4.5, 5.5, 1.5, 1.5), mex=0.5, cex.lab=1.5)	
	plot(type="n", xlim=xlim, ylim=ylim, x=NULL, y=NULL, xlab="V1", ylab="V2")

	# build grid data to predict
	datagrid <- gridGen(xlim, ylim, steps)
	
	nnet.pred <- predict(nnet.model, datagrid)
	# arg max
	nnet.pred <- max.col(nnet.pred)
	# avoid black
	nnet.pred <- nnet.pred + 1

	# set up grid viewport
	pos <- par("plt")
	grid::pushViewport(grid::viewport(x=grid::unit(pos[1], "npc"), y=grid::unit(pos[3], "npc"), 
		width=grid::unit(pos[2]-pos[1], "npc"), height=grid::unit(pos[4]-pos[3], "npc"), 
		just=c("left", "bottom")))
	grid::pushViewport(grid::viewport(xscale=xlim, yscale=ylim))


	# create alpha blended colors from labels
	#gridcolors <- numeric()
	#for(i in 1:dim(datagrid)[1]) {	
	#	curcolor <- getColor(nnet.pred[i])
	#	curcolor <- as.numeric(col2rgb(curcolor)) / 255
	#	curcolor <- rgb(curcolor[1], curcolor[2], curcolor[3], alpha)
	#	gridcolors[i] <- curcolor
	#}
	tempcolors <- getColor(nnet.pred)
	tempcolors <- col2rgb(tempcolors) / 255
	gridcolors <- rgb(tempcolors[1,], tempcolors[2,], tempcolors[3,], alpha)

	
	# plot decision regions
	grid::grid.rect(x=grid::unit(datagrid[,1], "native"), y=grid::unit(datagrid[,2], "native"), 
		width=grid::unit(xstep, "native"), height=grid::unit(ystep, "native"), 
		gp=grid::gpar(fill=gridcolors, lty=0))
		
	# then plot points, if not null
	if(displayPoints) {
		datalabels <- max.col(datalabels)
		labs <- as.numeric(datalabels[subset]) + 1
		newlabs <- character()
		for(i in 1:length(labs)) {
			newlabs[i] <- getColor(labs[i])
		}
		labs <- newlabs
		points(datamatrix[subset,1], datamatrix[subset,2], col=labs, pch=4, lwd=lwd)
	}	
	
	
}



