# $Id: NMBU.Graphics.R 37 2014-01-20 19:23:09Z khliland $

#################################
# Plotting for PCA/PCR/PLSR
principalComponentPlots <- function(){
    initializeDialog(title=gettextRcmdr("PCA/PCR/PLS plots"))
    variablesFrame <- tkframe(top)
	.numeric <- Variables()
	.activeModel <- ActiveModel()
	.activeDataSet <- ActiveDataSet()
    #.numeric <- Numeric()
	boxFrame1 <- tkframe(top)
	boxFrame2 <- tkframe(top)
	boxFrame3 <- tkframe(top)
	boxFrame4 <- tkframe(top)
    xBox <- variableListBox(boxFrame1, .numeric,
        title=gettextRcmdr("Colour/label points (select <= 1)"))
	compVar1 <- tclVar("1")
	compEntry1 <- ttkentry(variablesFrame, width="3", textvariable=compVar1)
	compVar2 <- tclVar("2")
	compEntry2 <- ttkentry(variablesFrame, width="3", textvariable=compVar2)

	model.class <- justDoIt(paste("class(", .activeModel, ")", sep=""))
	
    onOK <- function(){ # Actions to perform
		x <- getSelection(xBox)
        closeDialog()
        if (0 != length(x)) {
            # Ta vare p? og bruk til fargekoding
            }
		comp1 <- tclvalue(compVar1)
		comp2 <- tclvalue(compVar2)
		if(trim.blanks(comp1) == gettextRcmdr("") || trim.blanks(comp2) == gettextRcmdr("")){
			warning('Specify both components')
		}
		
        plottype <- as.character(tclvalue(plottypeVariable))
		the.axis.cross <- tclvalue(axis.crossVariable)
		if(plottype == gettextRcmdr("rmsep")){
			command <- paste("plot(RMSEP(", .activeModel, "))", sep="")
			logger(command)
			justDoIt(command)
		}
		else {
			if(plottype == gettextRcmdr("biplot")){
				if(model.class == "mvr"){
					command <- paste("biplot(", .activeModel, ", comps=c(", comp1, ",", comp2, "), main='Biplot')", sep="")
				} else {
					command <- paste("biplot(", .activeModel, ", choices=c(", comp1, ",", comp2, "), main='Biplot')", sep="")
				}
			}
			else {
				the.names <- tclvalue(namesVariable)
				if(plottype == gettextRcmdr("scores")){
					command <- paste("scoreplot(", .activeModel, ", main='Scoreplot', comps=c(", comp1, ",", comp2, ")", sep="")
					if(length(x) !=0){ # Colour coding by variable
						if(the.names == gettextRcmdr("1")){
							these.names <- paste(.activeDataSet, "[,'", x, "']", sep="") 
#							logger(paste("the.labels", " <- ", these.names, sep=""))
#							assign("the.labels", justDoIt(these.names), envir=.GlobalEnv)
							doItAndPrint(paste("the.labels", " <- ", these.names, sep=""))
							command <- paste(command, ", labels=the.labels, sub='Labels: ", x ,"'", sep="")
						}
						else{
							col.command <- paste("make.colours(",.activeDataSet, "[,'", x, "'])", sep="") 
#							logger(paste("the.colours", " <- ", col.command, sep=""))
#							assign("the.colours", justDoIt(col.command), envir=.GlobalEnv)
							doItAndPrint(paste("the.colours", " <- ", col.command, sep=""))
							command <- paste(command, ", col=the.colours, sub='Colours: ", x ,"'", sep="")
						}
					}
					else {
						if(the.names == gettextRcmdr("1"))
							command <- paste(command, ", labels='names', sub='Labels: observations'", sep="")
					}
				}
				else {
					the.spectra <- tclvalue(spectraVariable)
					if(plottype == gettextRcmdr("loadings")){
						if(the.spectra == gettextRcmdr("0"))
							command <- paste("loadingplot(", .activeModel, ", main='Loadingplot', comps=c(", comp1, ",", comp2, "), scatter=TRUE", sep="")
						else
							command <- paste("loadingplot(", .activeModel, ", main='Loadingplot', comps=c(", comp1, ",", comp2, "), scatter=FALSE", sep="")
					}
					if(plottype == gettextRcmdr("corr")){
						if(model.class == "mvr"){
							command <- paste("corrplot(", .activeModel, ", main='Correlation loadingplot', comps=c(", comp1, ",", comp2, ")", sep="")
						}
						else {
							doItAndPrint(paste("corr.load <- cor(",.activeDataSet,"[,rownames(loadings(", .activeModel, "))], scores(",.activeModel, "))",sep=""))
							#try(eval(parse(text=paste("n <- dim(loadings(",.activeModel,"))"))))
							#try(eval(parse(text=paste("X <- ",.activeDataSet,"[,rownames(loadings(", .activeModel, "))]", sep=""))))
							#CC <- matrix(NA,n[1],n[2])
							#for(i in 1:n[1]){ for(j in 1:n[2]){  try(eval(parse(text=paste("CC[i,j] <- cor(X[,i], scores(",.activeModel,")[,j])")))) }}
							#try(eval(parse(text=paste("colnames(CC) <- colnames(loadings(",.activeModel,"))"))))
							#try(eval(parse(text=paste("rownames(CC) <- rownames(loadings(",.activeModel,"))"))))
							#assign("corr.load", CC, envir=.GlobalEnv)
							command <- paste("corrplot(corr.load, main='Correlation loadingplot', comps=c(", comp1, ",", comp2, ")", sep="")
						}
					}
					if(the.names == gettextRcmdr("1"))
						command <- paste(command, ", labels='names', sub='Labels: variables'", sep="")
				}
				command <- paste(command, ")", sep="")
			}
			logger(command)
			justDoIt(command)
			if(the.axis.cross == gettextRcmdr("1")){
				abline(v=0, col="grey", lty=2)
				abline(h=0, col="grey", lty=2)
			}
		}
        tkfocus(CommanderWindow())
        }
	# Set up GUI
    OKCancelHelp(helpSubject="plot", model=TRUE)
	if(model.class == "mvr"){
		radioButtons(name="plottype", buttons=c("Loadings", "Scores", "Both", "Corr", "RMSEP"), values=c("loadings", "scores", "biplot", "corr", "rmsep"),
			labels=gettextRcmdr(c("Loadingplot", "Scoreplot", "Biplot", "Correlation loading plot", "RMSEP plot")), title=gettextRcmdr("Plot type"))}
	else{
		radioButtons(name="plottype", buttons=c("Loadings", "Scores", "Both", "Corr"), values=c("loadings", "scores", "biplot", "corr"),
			labels=gettextRcmdr(c("Loadingplot", "Scoreplot", "Biplot", "Correlation loading plot")), title=gettextRcmdr("Plot type"))}
 	tkgrid(plottypeFrame, row=1, column=1, sticky="w")
 	tkgrid(labelRcmdr(variablesFrame, text=gettextRcmdr("Component number for the x axis")), compEntry1, sticky="w")
	tkgrid(labelRcmdr(variablesFrame, text=gettextRcmdr("Component number for the y axis")), compEntry2, sticky="w")
 	tkgrid(variablesFrame, row=2, column=1, columnspan=2, sticky="w")
 	checkBoxes(frame="boxFrame2", boxes=c("names"), initialValues=c("0"),
        labels=gettextRcmdr(c("Use labels")))
 	checkBoxes(frame="boxFrame3", boxes=c("axis.cross"), initialValues=c("1"),
        labels=gettextRcmdr(c("Show axis cross")))
 	checkBoxes(frame="boxFrame4", boxes=c("spectra"), initialValues=c("0"),
        labels=gettextRcmdr(c("Spectra (loadings)")))
    tkgrid(getFrame(xBox), sticky="nw")
    tkgrid(boxFrame1, row=3, column=1, rowspan=3, sticky="w")
    tkgrid(boxFrame2, row=3, column=2, sticky="e")
    tkgrid(boxFrame3, row=4, column=2, sticky="e")
    tkgrid(boxFrame4, row=5, column=2, sticky="e")
    tkgrid(buttonsFrame, row=6, column=1, columnspan=2, stick="s")
    tkgrid.configure(helpButton, sticky="se")
    dialogSuffix(rows=6, columns=1)
}


#########################
## 2D discriminant plot
plotDA <- function(DAobject=NULL, regions=TRUE, contours=FALSE, resolution=100){
	if(is.null(DAobject)){
		try(eval(parse(text=paste("DAobject <- ", activeModel(), sep=""))))
	}

	data.name <- DAobject$call$data
	var.names <- colnames(DAobject$means) # Find and choose correct data
	X <- as.matrix(eval(parse(text = paste(data.name, "[,c('", paste(var.names, collapse="','"), "')]", sep=""))))
	n <- dim(X)
	if(n[2]>2)
		stop('Only two dimensions can be plotted (too many variables specified)')
	if(n[2]<2)
		stop('Only two dimensions can be plotted (too few variables specified)')
	# Initialisation
	y       <- eval(parse(text = paste(data.name, "[,c('", paste(formula(DAobject)[[2]], collapse="','"), "')]", sep="")))
	groups  <- levels(y)
	g       <- length(groups)
	n.g     <- numeric(g)
	index.g <- list()
	Y.dummy <- matrix(0,n[1],g)
	for(i in 1:g){ # Group-wise computations
		Y.dummy[,i]  <- y==groups[i]
		n.g[i]       <- sum(Y.dummy[,i])
		index.g[[i]] <- which(Y.dummy[,i]==1)
	}
	colnames(Y.dummy) <- groups
	names(n.g) <- groups

# Possible future dimension reduction for more than 2 columns of data.
#	## Reduce dimensions if >2
#	if(n[2] > 2){
#		Sy <- Y.dummy%*%inv(t(Y.dummy)%*%Y.dummy)%*%t(Y.dummy)
#		B  <- 1/n[1]*t(X)%*%Sy%*%X                      # Between groups covariance
#		To <- 1/(n[1]-1) * t(X)%*%(diag(rep(1,n[1]))-1/n[1]*matrix(1,n[1],n[1]))%*%X # Total centered covariance
#		#To <- 1/(n-1)*T'*T;                       # Total uncentered covariance
#		Wa <- To - B                               # Within groups covariance
#		Wa[abs(Wa)<.Machine$double.eps] <- 0
#		WB <- solve(Wa)%*%B
#		eigs <- eigen(WB)
#
#		XX <- X %*% eigs$vectors[,1]; XX <- XX%*%sign(XX[1])   # Discriminant vector 1 (training)
#		XY <- X %*% eigs$vectors[,2]; XY <- XY%*%sign(XY[1])   # Discriminant vector 2 (training)
#		X <- cbind(XX, XY)
#		colnames(X) <- c('BW1', 'BW2')
#	}

	#Estimating multinormal density parameters
	the.means <- matrix(0, min(n[2],2), g)
	the.covs  <- list()
	for(i in 1:g){ # Group-wise computations
		the.means[,i] <- apply(X[index.g[[i]],],2,mean)
		the.covs[[i]] <- var(X[index.g[[i]],])
	}
	S <- matrix(0, min(n[2],2),min(n[2],2))
	for(i in 1:g){ # Common covariance matrix
		S <- S + (n.g[i]-1)/(n[1]-g)*the.covs[[i]]
	}
	
	# Plotting limits
	min.x <- min(X[,1]); max.x <- max(X[,1]); dev.x <- max.x-min.x
	min.y <- min(X[,2]); max.y <- max(X[,2]); dev.y <- max.y-min.y
	y <- seq(min.y-0.15*dev.y,max.y+0.15*dev.y, length.out=resolution)
	x <- seq(min.x-0.15*dev.x,max.x+0.15*dev.x, length.out=resolution)
	nx <- length(x)
	ny <- length(y)

	# Matrices for holding probabilities over the grid spanned by x and y
	probs <- list()
	for(i in 1:g){
		probs[[i]] <- matrix(0,nx,ny)
	}
	ids <- matrix(0,nx,ny)

	# Computing multivariate density values over the grid spanned by x and y
	coords <- expand.grid(x,y)
	colnames(coords) <- colnames(X)
	preds <- predict(DAobject,coords)
	for(i in 1:g){
		probs[[i]] <- matrix(preds$posterior[,i],resolution,resolution)
	}
	if(regions==TRUE){ # Computations for region plotting
		for(i in 1:nx){
			for(j in 1:ny){
				id.temp <- numeric(g)
				for(k in 1:g){
					id.temp[k] <- probs[[k]][i,j]
				}
				ids[i,j] <- which(max(id.temp) == id.temp)
			}
		}
	}
	if(contours==TRUE){ # Computations for contour plotting
		if(class(DAobject)=='lda'){
			for(i in 1:g){
				probs[[i]] <- matrix(dmvnorm(coords,the.means[,i],S),resolution,resolution)
			}
		} else { # qda
			for(i in 1:g){
				probs[[i]] <- matrix(dmvnorm(coords,the.means[,i],the.covs[[i]]),resolution,resolution)
			}
		}
	} 

	# Using colors to indicate which density is the largest
	light.cols <- c('salmon', 'lightblue', 'lightgreen', 'lightgrey', 'yellow','tan','wheat')
	dark.cols <- c('red', 'blue', 'darkgreen', 'darkgrey', 'orange','brown','sienna')

	# Decision regions plot
	plot(median(x),median(y),type="n", xlim=c(min.x-0.1*dev.x,max.x+0.1*dev.x), ylim=c(min.y-0.1*dev.y,max.y+0.1*dev.y), xlab=colnames(X)[1], ylab=colnames(X)[2], main="Decision regions")
	if(regions==TRUE){
		colmat <- matrix("",nx,ny)
		for(i in 1:g){
			colmat[ids==i] <- light.cols[i]
		}
		for(i in 1:nx){
			for(j in 1:ny){
				points(x[i], y[j], pch=15, cex=2*40/(resolution), col=colmat[i,j] )
			}
		}
	}

	# Mean values for groups
	for(i in 1:g){
		points(the.means[2,i]~the.means[1,i],pch=17, cex=0.8, col=dark.cols[i])
		points(the.means[2,i]~the.means[1,i],pch=2, cex=0.8, col='black')
	}

	#Contour plots - densities
	if(contours==TRUE){
		for(i in 1:g){
			contour(x,y,probs[[i]], col=dark.cols[i], add=TRUE, drawlabels=FALSE, nlevels=7)
		}
	}

	# Observations
	for(i in 1:g){
		points(X[index.g[[i]],2]~X[index.g[[i]],1], col=dark.cols[i], pch=20, cex=2)
		points(X[index.g[[i]],2]~X[index.g[[i]],1], col='black', pch=21, cex=1.45) # Black filling
	}
}


############################
# dotPlot from package BHH2
dotPlot <- function (x, y = 0, xlim = range(x, na.rm = TRUE), xlab = NULL, 
    scatter = FALSE, hmax = 1, base = TRUE, axes = TRUE, frame = FALSE, 
    pch = 21, pch.size = "x", labels = NULL, hcex = 1, cex = par("cex"), 
    cex.axis = par("cex.axis"), ...) 
{
    if (is.null(xlab)) 
        xlab <- deparse(substitute(x))
    x <- x[!is.na(x)]
    xpd <- par("xpd")
    par(xpd = TRUE)
    on.exit(par(xpd = xpd))
    plot(c(0, 1), c(0, 1), xlim = xlim, type = "n", axes = FALSE, 
        cex = cex, cex.axis = cex.axis, frame = frame, xlab = xlab, 
        ylab = "", ...)
    if (axes) 
        axis(1, cex.axis = cex.axis)
    if (scatter) {
        dots(x, y = y, xlim = xlim, stacked = FALSE, hmax = hmax, 
            base = base, axes = FALSE, pch = pch, pch.size = pch.size, 
            labels = labels, hcex = hcex, cex = cex, cex.axis = cex.axis)
        y <- y + 2 * strheight(pch.size, units = "user")
        xlab <- ""
        axes <- FALSE
        base = FALSE
    }
    coord <- dots(x, y, xlim = xlim, stacked = TRUE, hmax = hmax, 
        base = base, axes = FALSE, pch = pch, pch.size = pch.size, 
        labels = labels, hcex = hcex, cex = cex, cex.axis = cex.axis)
    invisible(coord)
}
# ... and from BHH2
dots <- function (x, y = 0.1, xlim = range(x, na.rm = TRUE), stacked = FALSE, 
    hmax = 0.5, base = TRUE, axes = FALSE, pch = 21, pch.size = "x", 
    labels = NULL, hcex = 1, cex = par("cex"), cex.axis = par("cex.axis")) 
{
    x <- x[!is.na(x)]
    hdots <- y
    xmin <- xlim[1]
    xmax <- xlim[2]
    x <- x[(x >= xmin) & (x <= xmax)]
    b <- strwidth(pch.size, units = "user", cex = cex)
    h <- strheight(pch.size, units = "user", cex = hcex * cex)
    if (stacked) {
        if (xmax - xmin < b) {
            stop("x-dimension resolution problem")
        }
        else {
            xu <- seq(xmin, xmax, by = b)
        }
        m <- length(xu)
        tab <- data.frame(j = 1:m, k = rep(0, m), xu = xu)
        n <- length(x)
        y <- rep(0, n)
        for (i in 1:n) {
            l <- max(tab$j[tab$xu <= x[i]])
            x[i] <- xu[l] + b/2
            tab$k[l] <- 1 + tab$k[l]
            y[i] <- tab$k[l]
        }
        y <- y * h
        u <- hdots + max(y)
        if (hmax <= hdots) 
            warning(paste("dot base <hdots=", hdots, "> higher than maximum column height <hmax=", 
                hmax, ">...", sep = ""))
        if (u > hmax) 
            y <- (hmax - hdots) * y/u
        y <- hdots + y
    }
    else {
        y <- rep(hdots, length(x))
    }
    if (!is.null(labels)) 
        text(x, y, labels = labels, cex = cex)
    else points(x, y, pch = pch, cex = cex)
    points.coord <- data.frame(x, y)
    if (axes) {
        segments(xmin - b, hdots - h/4, xmax + b, hdots - h/4)
        x <- pretty(x, n = 3, h = 0.5)
        x <- x[(x > xmin - b) & (x < xmax + b)]
        for (i in seq(x)) segments(x[i], hdots - h/4, x[i], hdots - 
            h/2)
        y <- rep(hdots - h, length(x))
        text(x, y, labels = x, cex = cex.axis)
    }
    if (base && !axes) 
        segments(xmin - b, hdots - h/4, xmax + b, hdots - h/4)
    invisible(points.coord)
}


############################
# Contour plot of mixture designs (three factors/components)
mixture.contour <- function(data, formula1, n.tick=6, n.grade=15, resolution=3, FUN=NULL, FUN2=NULL, PLS=FALSE, ncomp=NULL,
	mix.format="dec", show.points=FALSE, show.contour=TRUE, zoomed=FALSE, pch=21, cex=1.25, col.points="black", fill.points="white"){
	if(show.contour){
		data <- data[,c(attr(terms(formula1),"term.labels")[1:3],setdiff(all.vars(formula1),attr(terms(formula1),"term.labels")[1:3]))]
	} else {
		data <- data[,c(attr(terms(formula1),"term.labels")[1:3])]
		data$resp <- rnorm(dim(data)[1])
		formula1 <- update(formula1,resp~.)
	}
	data.names <- colnames(data)

	orders <- 1:3

	is100 <- TRUE
	if(round(sum(data[1,1:3]))==1){
#		data[,1:3] <- data[,1:3]*100
		is100 <- FALSE
	}
	if(zoomed){
		col.min <- apply(data[,-4],2,min)
		pseudo  <- t(t(data[,-4])-col.min)
		row.sum <- sum(pseudo[1,])
		pseudo  <- pseudo/row.sum
		col.max <- row.sum+col.min
	} else {
		col.min <- rep(0,3)
		if(is100){
			pseudo  <- data[,-4]/100
			col.max <- rep(100,3)
		} else {
			col.max <- rep(1,3)
			pseudo  <- data[,-4]
		}
	}
	
	pseudo <- as.data.frame(cbind(pseudo,data[,4]))
	colnames(pseudo) <- c(data.names[orders],data.names[4])
	
	if(PLS){
		m1 <- plsr(formula1, data=pseudo, ncomp=ncomp)
	} else {
		m1 <- lm(formula1, data=pseudo)}
	
	if(!show.contour)
		resolution <- 1
	trian <- expand.grid(base=seq(0,1,l=100*resolution), high=seq(0,sin(pi/3),l=87*resolution))
	trian <- subset(trian, (base*sin(pi/3)*2)>high)
	trian <- subset(trian, ((1-base)*sin(pi/3)*2)>high)

	new2 <- data.frame(a=trian$high*2/sqrt(3))
	new2$b <- trian$base-trian$high/sqrt(3)
	new2$c <- 1-new2$b-new2$a
	colnames(new2) <- data.names[orders[c(2,3,1)]]
	
	if(PLS){
		trian$yhat <- predict(m1, newdata=new2)[,1,ncomp]
	} else {
		trian$yhat <- predict(m1, newdata=new2)
	}
	if(show.contour){
		Gray <- function(n){gray(seq(0,1, length.out=n))}
	} else {
		Gray <- function(n){gray(rep(1,n))}
	}
	my.frac <- function(inn){
		hel <- floor(inn)
		dec <- inn-hel
		out <- as.character(hel)
		for(i in 1:length(out)){
			if(dec[i]>0){
				if(hel[i]>0){
					out[i] <- paste(out[i], fractions(dec[i]))
				} else {
					out[i] <- paste(fractions(dec[i]))
				}
			}
		}
		out
	}
	format.mix <- function(mix.format, inn) {
		switch(mix.format,
		dec = inn,
		frac = my.frac(inn),
		ratio = round(inn*(n.tick-1))
		)}
		
	grade.trellis <- function(pseudo, data, col.min, col.max, n.tick, col=1, lty=2, lwd=0.8, col.points=col.points, fill.points=fill.points, mix.format, show.points, cex, pch){
	  tt <- rbind(col.min,col.max)
	  t1 <- tt[,1]; t2 <- tt[,2]; t3 <- tt[,3]
	  t1 <- seq(t1[1],t1[2],length.out=n.tick)[2:(n.tick-1)]; t2 <- seq(t2[1],t2[2],length.out=n.tick)[2:(n.tick-1)]; t3 <- seq(t3[1],t3[2],length.out=n.tick)[2:(n.tick-1)]
	  t1 <- format.mix(mix.format,t1); t2 <- format.mix(mix.format,t2); t3 <- format.mix(mix.format,t3)
	  x1 <- seq(0, 1, length.out=n.tick)[2:(n.tick-1)]
	  x2 <- x1/2
	  y2 <- x1*sqrt(3)/2
	  x3 <- (1-x1)*0.5+x1
	  y3 <- sqrt(3)/2-x1*sqrt(3)/2
	  panel.segments(x1, 0, x2, y2, col=col, lty=lty, lwd=lwd)
	  panel.text(x1, 0, label=t3, pos=1)
	  panel.segments(x1, 0, x3, y3, col=col, lty=lty, lwd=lwd)
	  panel.text(x2, y2, label=rev(t1), pos=2)
	  panel.segments(x2, y2, 1-x2, y2, col=col, lty=lty, lwd=lwd)
	  panel.text(x3, y3, label=rev(t2), pos=4)
	  if(show.points){
  		panel.points(1-pseudo[,1]-pseudo[,2]/2,pseudo[,2]*sqrt(3)/2,col=col.points, fill=fill.points, pch=pch, cex=cex, lwd=2)
	  }
	}

	if(is.null(FUN)){
		a <- b <- pretty(seq(min(trian$yhat), max(trian$yhat), length.out=n.grade), n=n.grade)# min.n=n.grade %/% 1.5)
	} else {
		b <- 10^pretty(log10(lapply(list(seq(min(trian$yhat), max(trian$yhat), length.out=n.grade)),FUN)[[1]]), n=n.grade)
		a <- lapply(list(b),FUN2)[[1]]
	}
	trihat.values <- list(a=a,b=b)
	n.grade <- length(trihat.values$a)-1

	if(show.contour){
		print(levelplot(yhat~base*high, trian, aspect="iso", xlim=c(-0.1,1.1), ylim=c(-0.1,0.96), cuts=n.grade-1,
				  xlab=NULL, ylab=NULL, contour=TRUE, col.regions = Gray, colorkey=list(at=trihat.values$a,labels=format(trihat.values$b,digits=1,nsmall=1,scientific=FALSE), col=Gray(n.grade), axis.line=list(col=1)),
				  par.settings=list(axis.line=list(col=NA), axis.text=list(col=1)), scales = list(draw=FALSE),
				  panel=function(..., at, contour=TRUE, labels=NULL){
					panel.levelplot(..., at=at, contour=contour, 
									lty=3, lwd=0.5, col=1)
				  }))
	} else {
		print(levelplot(yhat~base*high, trian, aspect="iso", xlim=c(-0.1,1.1), ylim=c(-0.1,0.96), cuts=n.grade-1,
				  xlab=NULL, ylab=NULL, contour=TRUE, col.regions = Gray, colorkey=NULL,
				  par.settings=list(axis.line=list(col=NA), axis.text=list(col=1)), scales = list(draw=FALSE),
				  panel=function(..., at, contour=TRUE, labels=NULL){
					panel.levelplot(..., at=at, contour=contour, 
									lty=3, lwd=0.5, col="white")
				  }))
	}
	trellis.par.set(list(clip = list(panel=FALSE)))
	trellis.focus("panel", 1, 1, highlight=FALSE)
	panel.segments(c(0,0,0.5), c(0,0,sqrt(3)/2), c(1,1/2,1), c(0,sqrt(3)/2,0), lwd=2)
	grade.trellis(pseudo, data, col.min, col.max, n.tick, mix.format=mix.format, show.points=show.points, cex=cex, pch=pch, col.points=col.points, fill.points=fill.points)
	panel.text(0, 0, label=data.names[orders[1]], pos=2)
	panel.text(1/2, sqrt(3)/2, label=data.names[orders[2]], pos=3)
	panel.text(1, 0, label=data.names[orders[3]], pos=4)
	trellis.unfocus()
}

#####################################
# Fitted line plot
CIplot <- function (lm.object, xlim = range(data[, x.name]), newdata, conf.level = 0.95, 
    data = model.frame(lm.object), newfit, ylim, pch = 16, main.cex = 1, 
    main = list(paste(100 * conf.level, "% confidence and prediction intervals, R^2=", format(summary(lm.object)[]$r.squared, digits=2), sep = ""), cex = main.cex), ...) 
{
	lm.coef <- coef(lm.object)
	if(length(lm.coef)>1){
		fit.text <- paste(format(lm.coef[1],digits=2), " + ", format(lm.coef[2],digits=2), "*", names(lm.coef)[2], sep="")}
	else{
		fit.text <- paste(format(lm.coef[1],digits=2), "*", names(lm.coef)[1], sep="")}		

    formula.lm <- formula(lm.object)
    x.name <- as.character(formula.lm[[3]])
    y.name <- as.character(formula.lm[[2]])
    missing.xlim <- missing(xlim)
    missing.ylim <- missing(ylim)
    missing.newdata <- missing(newdata)
    # if.R(s = {
        # my.data.name <- as.character(lm.object$call$data)
        # if (length(my.data.name) == 0) 
            # stop("Please provide an lm.object calculated with an explicit 'data=my.data.frame' argument.")
        # undo.it <- (!is.na(match(my.data.name, objects(0))))
        # if (undo.it) 
            # old.contents <- get(my.data.name, frame = 0)
        # my.data <- try(get(my.data.name))
        # if (class(my.data) == "Error") 
            # my.data <- try(get(my.data.name, frame = sys.parent()))
        # if (class(my.data) == "Error") 
            # stop("Please send me an email with a reproducible situation that got you here. (rmh@temple.edu)")
        # assign(my.data.name, my.data, frame = 0)
    # }, r = {
    # })
    default.newdata <- data.frame(seq(xlim[1], xlim[2], length = 51))
    names(default.newdata) <- x.name
    if (missing.xlim) 
        xlim <- xlim + diff(xlim) * c(-0.02, 0.02)
    if (missing.newdata) {
        newdata <- default.newdata
        newdata.x <- numeric()
    }
    else {
        if (is.na(match(x.name, names(newdata)))) 
            stop(paste("'newdata' must be a data.frame containing a column named '", 
                x.name, "'", sep = ""))
        if (missing.xlim) 
            xlim = range(xlim, newdata[[x.name]])
        newdata.x <- as.data.frame(newdata)[, x.name]
        newdata <- rbind(as.data.frame(newdata)[, x.name, drop = FALSE], 
            default.newdata)
        newdata <- newdata[order(newdata[, x.name]), , drop = FALSE]
    }
    if (missing.xlim) 
        xlim <- xlim + diff(xlim) * c(-0.02, 0.02)
    if (missing(newfit)) 
        newfitF <- function(){
            new.p <- predict(lm.object, newdata = newdata, se.fit = TRUE, 
                level = conf.level, interval = "prediction")
            new.c <- predict(lm.object, newdata = newdata, se.fit = TRUE, 
                level = conf.level, interval = "confidence")
            tmp <- new.p
            tmp$ci.fit <- new.c$fit[, c("lwr", "upr"), drop = FALSE]
            dimnames(tmp$ci.fit)[[2]] <- c("lower", "upper")
            attr(tmp$ci.fit, "conf.level") <- conf.level
            tmp$pi.fit <- new.p$fit[, c("lwr", "upr"), drop = FALSE]
            dimnames(tmp$pi.fit)[[2]] <- c("lower", "upper")
            attr(tmp$pi.fit, "conf.level") <- conf.level
            tmp$fit <- tmp$fit[, "fit", drop = FALSE]
            tmp
        }
	newfit <- newfitF()
    tpgsl <- trellis.par.get("superpose.line")
    tpgsl <- Rows(tpgsl, 1:4)
    tpgsl$col[1] <- 0
    if (missing.ylim) {
        ylim <- range(newfit$pi.fit, data[, y.name])
        ylim <- ylim + diff(ylim) * c(-0.02, 0.02)
    }
    xyplot(formula.lm, data = data, newdata = newdata, newfit = newfit, 
        newdata.x = newdata.x, xlim = xlim, ylim = ylim, pch = pch, 
        panel = function(..., newdata.x) {
            panel.ci.plot(...)
            if (length(newdata.x) > 0) 
                panel.rug(x = newdata.x)
        }, main = main, key = list(border = TRUE, space = "right", 
            text = list(c("observed", fit.text, "conf int", "pred int")), 
            points = list(pch = c(pch, 32, 32, 32), col = c(trellis.par.get("plot.symbol")$col, 
                tpgsl$col[2:4])), lines = tpgsl), ...)
}
 panel.ci.plot <- function (x, y, newdata, newfit = newfit, ...) 
{
    tpgsl <- trellis.par.get("superpose.line")
    panel.xyplot(x, y, ...)
    panel.xyplot(x = newdata[[1]], y = newfit$fit, type = "l", 
        lty = tpgsl$lty[2], col = tpgsl$col[2], lwd = tpgsl$lwd[2])
    panel.xyplot(x = newdata[[1]], y = newfit$ci.fit[, "lower"], 
        type = "l", lty = tpgsl$lty[3], col = tpgsl$col[3], lwd = tpgsl$lwd[3])
    panel.xyplot(x = newdata[[1]], y = newfit$ci.fit[, "upper"], 
        type = "l", lty = tpgsl$lty[3], col = tpgsl$col[3], lwd = tpgsl$lwd[3])
    panel.xyplot(x = newdata[[1]], y = newfit$pi.fit[, "lower"], 
        type = "l", lty = tpgsl$lty[4], col = tpgsl$col[4], lwd = tpgsl$lwd[4])
    panel.xyplot(x = newdata[[1]], y = newfit$pi.fit[, "upper"], 
        type = "l", lty = tpgsl$lty[4], col = tpgsl$col[4], lwd = tpgsl$lwd[4])
    if.R(s = {
        axis(1, at = mean(x), tck = -0.03, labels = FALSE)
        axis(1, at = mean(x), ticks = FALSE, labels = "xbar", 
            line = 0.9)
        axis(3, at = mean(x), tck = -0.03, labels = FALSE)
        axis(3, at = mean(x), ticks = FALSE, labels = "xbar")
    }, r = {
        cpl <- current.panel.limits()
        pushViewport(viewport(xscale = cpl$xlim, yscale = cpl$ylim, 
            clip = "off"))
        panel.axis("bottom", at = mean(x), tck = 1.5, labels = FALSE)
        panel.axis("bottom", at = mean(x), ticks = FALSE, labels = expression(bar(x)), 
            rot = 0, outside = TRUE)
        panel.axis("top", at = mean(x), tck = 1.5, labels = FALSE)
        panel.axis("top", at = mean(x), ticks = FALSE, labels = expression(bar(x)), 
            rot = 0, outside = TRUE)
        upViewport()
    })
}
if.R <- function (r, s) 
{
    if (exists("is.R") && is.function(is.R) && is.R()) 
        r
    else s
}


#####################################
# Customized basic diagnostics plot
plotModelNMBU <- function(){
    .activeModel <- ActiveModel()
    if (is.null(.activeModel) || !checkMethod("plot", .activeModel)) return()
    command <- "oldpar <- par(oma=c(0,0,3,0), mfrow=c(2,2))"
    justDoIt(command)
    logger(command)
    doItAndPrint(paste("plot(", .activeModel, ", which = c(1,2,4,5), add.smooth=FALSE)", sep=""))
    command <- "par(oldpar)"
    justDoIt(command)
    logger(command)
    }

#####################################
# Plot by groups
plotByGroups <- function(x.name, y.name, z.name, lineType, axisLabel, legend, col, cex=1, ...){
	.tmp.by.groups <- eval(parse(text=ActiveDataSet()))
	xLabel <- x.name
#	attach(.tmp.by.groups,name="tmp")
#	on.exit(detach("tmp"))
	n <- dim(.tmp.by.groups)[1]
	levs <- as.character(eval(parse(text=paste("unique(", ActiveDataSet(), '$', z.name, ")", sep=""))))
	if(missing(col))
		col <- 1:length(levs)
	if(length(col)==1)
		col <- rep(col,length(levs))
	xRange <- eval(parse(text=paste("range(", ActiveDataSet(), '$', x.name, ")", sep="")))
	yRange <- eval(parse(text=paste("range(", ActiveDataSet(), '$', y.name, ")", sep="")))
    top <- 3.5 + length(levs)
	if(legend){
		mar <- par("mar")
		eval(parse(text=paste(".mar <- par(mar=c(", mar[1], ",", mar[2], ",", top, ",", mar[4], "))", sep="")))
	}
	if(lineType=="p"){
		x1 <- eval(parse(text=paste("order(", ActiveDataSet(), '$', x.name, "[", ActiveDataSet(), '$', z.name, "=='", levs[1], "'])", sep="")))
		eval(parse(text=paste("plot(", ActiveDataSet(), '$', x.name, "[", ActiveDataSet(), '$', z.name, "=='", levs[1], "'][x1]", ", ", ActiveDataSet(), '$', y.name, "[", ActiveDataSet(), '$', z.name, "=='", levs[1], "'][x1]", ", type='p', col=col[1], cex=cex, xlab='", xLabel, "', ylab='", axisLabel, "', xlim=c(",xRange[1], ",", xRange[2],"), ylim=c(",yRange[1], ",", yRange[2],"), ...)", sep="")))
		for(i in 2:length(levs)){
			x1 <- eval(parse(text=paste("order(", ActiveDataSet(), '$', x.name, "[", ActiveDataSet(), '$', z.name, "=='", levs[i], "'])", sep="")))
			eval(parse(text=paste("points(", ActiveDataSet(), '$', x.name, "[", ActiveDataSet(), '$', z.name, "=='", levs[i], "'][x1]", ", ", ActiveDataSet(), '$', y.name, "[", ActiveDataSet(), '$', z.name, "=='", levs[i], "'][x1], col=col[", i,"], cex=cex, pch=",i,")", sep="")))
		}
	}
	if(lineType=="l" || lineType=="b"){
		x1 <- eval(parse(text=paste("order(", ActiveDataSet(), '$', x.name, "[", ActiveDataSet(), '$', z.name, "=='", levs[1], "'])", sep="")))
		eval(parse(text=paste("plot(", ActiveDataSet(), '$', x.name, "[", ActiveDataSet(), '$', z.name, "=='", levs[1], "'][x1]", ", ", ActiveDataSet(), '$', y.name, "[", ActiveDataSet(), '$', z.name, "=='", levs[1], "'][x1]", ", type='l', col=col[1], cex=cex, xlab='", xLabel, "', ylab='", axisLabel, "', xlim=c(",xRange[1], ",", xRange[2],"), ylim=c(",yRange[1], ",", yRange[2],"), ...)", sep="")))
		for(i in 2:length(levs)){
			x1 <- eval(parse(text=paste("order(", ActiveDataSet(), '$', x.name, "[", ActiveDataSet(), '$', z.name, "=='", levs[i], "'])", sep="")))
			eval(parse(text=paste("lines(", ActiveDataSet(), '$', x.name, "[", ActiveDataSet(), '$', z.name, "=='", levs[i], "'][x1]", ", ", ActiveDataSet(), '$', y.name, "[", ActiveDataSet(), '$', z.name, "=='", levs[i], "'][x1], col=col[", i,"], cex=cex,)", sep="")))
		}
	}
	if(lineType=="b"){
		for(i in 1:length(levs)){
			x1 <- eval(parse(text=paste("order(", ActiveDataSet(), '$', x.name, "[", ActiveDataSet(), '$', z.name, "=='", levs[i], "'])", sep="")))
			eval(parse(text=paste("points(", ActiveDataSet(), '$', x.name, "[", ActiveDataSet(), '$', z.name, "=='", levs[i], "'][x1]", ", ", ActiveDataSet(), '$', y.name, "[", ActiveDataSet(), '$', z.name, "=='", levs[i], "'][x1], col=col[", i,"], cex=cex,, pch=",i,")", sep="")))
		}
	}
	if(legend){
        .xpd <- par(xpd=TRUE)
        usr <- par("usr")
		eval(parse(text=paste("legend(", usr[1], ", ", usr[4] + 1.2*top*strheight("x.name"), ", legend=",
                paste("c(", paste(paste('"', levs, '"', sep=""), collapse=","), ")", sep=""),
                ", col=c(", paste(col, collapse=","), ")", ifelse(lineType=="l"||lineType=="b",", lty=1",""), ifelse(lineType=="p"||lineType=="b",paste(", pch=1:",length(levs),sep=""),""),")", sep="")))
        par(xpd=.xpd)
		par(mar=mar)
	}
}

#####################################
# Pie chart (by summaries)
pieChartNMBU <- function(){
	initializeDialog(title=gettextRcmdr("Pie chart from summaries"))
	labelsBox <- variableListBox(top, Factors(), title=gettextRcmdr("Labels (pick one)"))
	sizeBox <- variableListBox(top, Numeric(), title=gettextRcmdr("Sizes (pick one)"))
		optionsFrame <- tkframe(top)
#		pairwiseVariable <- tclVar("0")
#		pairwiseCheckBox <- tkcheckbutton(optionsFrame, variable=pairwiseVariable)
	onOK <- function(){
		labels <- getSelection(labelsBox)
		sizes <- getSelection(sizeBox)
		plottype <- as.character(tclvalue(plottypeVariable))
		closeDialog()
		if (length(labels) == 0){
			errorCondition(recall=pieChartNMBU, message=gettextRcmdr("You must select a variable for labels."))
			return()
		}
		if (length(sizes) == 0){
			errorCondition(recall=pieChartNMBU, message=gettextRcmdr("You must select a variable for sizes."))
			return()
		}

		.activeDataSet <- ActiveDataSet()
		n <- justDoIt(paste("length(", .activeDataSet, "$", labels, ")"))
		if(plottype == gettextRcmdr("pie")){
			command <- paste("pie(", .activeDataSet, "$", sizes, ", ", .activeDataSet, "$", labels, ", col=rainbow(", n,"))", sep="")
		} else {
			command <- paste("barplot(", .activeDataSet, "$", sizes, ", names.arg=", .activeDataSet, "$", labels, ", col=rainbow(", n,"))", sep="")
		}
		justDoIt(command)
		logger(command)
		tkfocus(CommanderWindow())
	}
	OKCancelHelp(helpSubject="pie", model=TRUE)
	tkgrid(getFrame(labelsBox), getFrame(sizeBox), sticky="nw")
#		tkgrid(labelRcmdr(optionsFrame, text=gettextRcmdr("Pairwise comparisons of means")), pairwiseCheckBox, sticky="w")
#		tkgrid(optionsFrame, sticky="w", columnspan=2)
	radioButtons(name="plottype", buttons=c("pie", "pie"), values=c("pie", "bar"),
		labels=gettextRcmdr(c("Pie chart", "Bar chart")), title=gettextRcmdr("Plot type"))
	tkgrid(plottypeFrame, columnspan=2, sticky="w")
	tkgrid(buttonsFrame, columnspan=2, sticky="w")
	dialogSuffix(rows=3, columns=2)
}


#####################################
# Add fitted density to plot (lines)
plotFitDens <- function(object, distr, scaling=1, col=2, lwd=2, ...){
	fits <- fitdistr(object,distr)$estimate
	varRange <- range(object)
	varRange <- c(varRange[1] - 0.5*diff(varRange), varRange[2] + 0.5*diff(varRange))
	x <- seq(varRange[1],varRange[2],length.out=100)
	if(distr == "normal"){
		lines(x, dnorm(x, fits[1], fits[2])*scaling, col=col, lwd=lwd, ...)
	}
	if(distr == "t"){
		lines(x, 1/fits[2] * dt((x-fits[1])/fits[2], fits[3])*scaling, col=col, lwd=lwd, ...)
	}
	if(distr == "exponential"){
		lines(x, dexp(x, fits[1])*scaling, col=col, lwd=lwd, ...)
	}
	if(distr == "gamma"){
		lines(x, dgamma(x, fits[1], fits[2])*scaling, col=col, lwd=lwd, ...)
	}
	if(distr == "geometric"){
		lines(round(x), dgeom(round(x), fits[1])*scaling, col=col, lwd=lwd, ...)
	}
	if(distr == "log-normal" || distr == "lognormal"){
		lines(x, dlnorm(x, fits[1], fits[2])*scaling, col=col, lwd=lwd, ...)
	}
	if(distr == "logistic"){
		lines(x, dlogis(x, fits[1], fits[2])*scaling, col=col, lwd=lwd, ...)
	}
	if(distr == "negative binomial"){
		lines(round(x), dnbinom(round(x), fits[1], mu=fits[2])*scaling, col=col, lwd=lwd, ...)
	}
	if(distr == "Poisson"){
		lines(round(x), dpois(round(x), fits[1])*scaling, col=col, lwd=lwd, ...)
	}
	if(distr == "weibull"){
		lines(x, dweibull(x, fits[1], fits[2])*scaling, col=col, lwd=lwd, ...)
	}
}


#####################################
# Discrete distribution plot with histograms
binomialDistributionPlotNMBU <- function(){discreteDistributionPlotNMBU("binomial")}
PoissonDistributionPlotNMBU <- function(){discreteDistributionPlotNMBU("Poisson")}
geomDistributionPlotNMBU  <- function(){discreteDistributionPlotNMBU("geom")}
hyperDistributionPlotNMBU  <- function(){discreteDistributionPlotNMBU("hyper")}
negbinomialDistributionPlotNMBU  <- function(){discreteDistributionPlotNMBU("negbinomial")}
discreteDistributionPlotNMBU <- function(nameVar){
#	fVar<-get(paste(nameVar,"Distribution",sep=""))
	fVar <- eval(parse(text=paste("Rcmdr:::",nameVar,"Distribution",sep="")))
	nnVar<-length(fVar$params)
	dialogName <- paste(nameVar,"DistributionPlot", sep="")
	defaults <- list(initialValues=fVar$initialValues, type="Probability")
	initial <- getDialog(dialogName, defaults=defaults)
	initializeDialog(title=gettextRcmdr(paste(fVar$titleName,"Distribution",sep=" ")))
	paramsVar<-paste(fVar$params,"Var",sep="")
	paramsEntry<-paste(fVar$params,"Entry",sep="")
	for (i in 1:nnVar) {
		eval(parse(text=paste(paramsVar[i],"<-tclVar('",initial$initialValues[i],"')",sep="")))
		eval(parse(text=paste(paramsEntry[i],"<-ttkentry(top, width='6', textvariable=",paramsVar[i],")",sep="")))
	}
	functionVar <- tclVar(initial$type)
	densityButton <- ttkradiobutton(top, variable=functionVar, value="Probability")
	distributionButton <- ttkradiobutton(top, variable=functionVar, value="Cumulative Probability")
	onOK <- function(){
#		nameVarF<-get(paste(nameVar,"DistributionPlot",sep=""),mode="function")
		nameVarF <- eval(parse(text=paste("Rcmdr:::",nameVar,"DistributionPlot",sep="")))
		closeDialog()
		warn <- options(warn=-1)
		vars <- double(nnVar)
		for (i in 1:nnVar) {
			vars[i]<-as.numeric(tclvalue(get(paramsVar[i])))
		}
		if (length(fVar$paramsRound)>0) {
			for (j in fVar$paramsRound) {
				vars[j]<-round(vars[j])
			}
		}
		options(warn)
		for (i in 1:length(fVar$errorConds)) {
			if (eval(parse(text=fVar$errorConds[i]))) {
				errorCondition(recall=nameVarF, message=gettextRcmdr(fVar$errorTexts[i]))
				return()
			}
		}
		fun <- tclvalue(functionVar)
		pasteVar<-""
		for (i in 1:nnVar) {
			pasteVar<-paste(pasteVar,", ",fVar$params[i],"=",vars[i],sep="")
		}
		min <- eval(parse(text=paste("q",fVar$funName,"(.0005",pasteVar,")",sep="")))
		max <- eval(parse(text=paste("q",fVar$funName,"(.9995",pasteVar,")",sep="")))
		command <- paste(min, ":", max, sep="")
#		logger(paste(".x <- ", command, sep=""))
#		assign(".x", justDoIt(command), envir=.GlobalEnv)
		doItAndPrint(paste(".x <- ", command, sep=""))
		switch(nameVar,
				"binomial" = xlabVar<-"Number of Successes",
				"Poisson" = xlabVar<-"x",
				"geom" = xlabVar<-"Number of Failures until Success",
				"hyper" = xlabVar<-"Number of White Balls in Sample",
				"negbinomial" = xlabVar <-"Number of Failures Until Target Successes"
		)
		mainVar<-""
		if (nameVar=="negbinomial") {
			mainVar<-paste(", Trials=",vars[1],", Prob=",vars[2],sep="")
		} else if (nameVar=="hyper") {
			mainVar<-paste(", m=",vars[1],", n=",vars[2],", k=",vars[3],sep="")
		} else {
			for (i in 1:nnVar) {
				mainVar<-paste(mainVar,", ", fVar$paramsLabels[i],"=",vars[i],sep="")
			}   
		}
		if (fun == "Probability"){
			doItAndPrint(paste("barplot(d",fVar$funName,"(.x", pasteVar,
							'), names.arg=.x, space=0, xlab="',xlabVar,'", ylab="Probability Mass", main="',fVar$titleName,
							' Distribution: ',substr(mainVar,2,nchar(mainVar)),'")', sep=""))
#			doItAndPrint(paste("points(.x, d",fVar$funName,"(.x", pasteVar,
#							'), pch=16)', sep=""))
		}
		else {
			command <- "rep(.x, rep(2, length(.x)))"
#			logger(paste(".x <- ", command, sep=""))
#			assign(".x", justDoIt(command), envir=.GlobalEnv)
			doItAndPrint(paste(".x <- ", command, sep=""))
			doItAndPrint(paste("plot(.x[-1], p",fVar$funName,"(.x",
							pasteVar,')[-length(.x)], xlab="',xlabVar,
							'",ylab="Cumulative Probability", main="',
							fVar$titleName,' Distribution: ',substr(mainVar,2,nchar(mainVar)),'", type="l")', sep=""))
		}
		doItAndPrint('abline(h=0, col="gray")')
		remove(.x, envir=.GlobalEnv)
		logger("remove(.x)")
		tkfocus(CommanderWindow())
		putDialog(dialogName, list(initialValues=vars, type=fun))#, resettable=FALSE)
	}
	OKCancelHelp(helpSubject=paste("d",fVar$funName,sep=""), reset=dialogName)
	for (i in 1:nnVar) {
		tkgrid(labelRcmdr(top, text=gettextRcmdr(fVar$paramsLabels[i])), get(paramsEntry[i]), sticky="e")
	}
	tkgrid(labelRcmdr(top, text=gettextRcmdr("Plot probability mass function")), densityButton, sticky="e")
	tkgrid(labelRcmdr(top, text=gettextRcmdr("Plot distribution function")), distributionButton, sticky="e")
	tkgrid(buttonsFrame, columnspan=2, sticky="w")
	for (i in 1:nnVar) {
		tkgrid.configure(get(paramsEntry[i]), sticky="w")
	}
	tkgrid.configure(densityButton, sticky="w")
	tkgrid.configure(distributionButton, sticky="w")
	dialogSuffix(rows=5, columns=2, focus=get(paramsEntry[1]))
}
