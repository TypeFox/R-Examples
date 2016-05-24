#===========================================================
# Plotting functions (in alphabetic order)                 |
#  addArrows    : Add arrows to current plot
#  addLabel     : Add label to current plot
#  addLegend    : Add a legend to current plot
#  drawBars     : Draw a linear barplot on the current graph
#  expandGraph  : Tweaks values to expand margins for multiple graphs
#  lucent       : Create translucent colour
#  pickCol      : Display interactive colour picking palette
#  plotACF      : Plot auto correlations (support for BRugs)
#  plotAsp      : Plots x and y vectors with plot() but maintaining a fixed aspect
#  plotBubbles  : Function to construct a bubble plot for a matrix z
#  plotCsum     : Plots cumulative frequecy of data
#  plotDens     : Plot density curves (support for BRugs)
#  plotFriedEggs: Pairs plot featuring fried eggs and beer
#  plotSidebars : Plot (x,y) table as horizontal sidebars
#  plotTrace    : Plot trace lines (support for BRugs)
#  resetGraph   : Resets par() values to R default
#  testAlpha    : Display various alpha transparencies
#  testCol      : Display test colours as circular patches
#  testLty      : Display line types available
#  testLwd      : Display line widths
#  testPch      : Display plotting symbols or octal strings
#-----------------------------------------------------------
# Authors:                                                 |
#  Jon T. Schnute <Jon.Schnute@dfo-mpo.gc.ca>,             |
#  Alex Couture-Beil <alex@mofo.ca>,                       |
#  Rowan Haigh <Rowan.Haigh@dfo-mpo.gc.ca>                 |
#  Rob Kronlund <Allen.Kronlund@dfo-mpo.gc.ca>             |
#===========================================================


#addArrows------------------------------2007-08-22
# Calls 'arrows' function using relative (0:1) coordinates
# Arguments:
#  x1 - draw from
#  y1 - draw from
#  x2 - draw to
#  y2 - draw to
#  ... - arguments used by key, such as "lines", "text", or "rectangle"
#-------------------------------------------JTS/RH
addArrows <- function (x1, y1, x2, y2, ...) {
	uxy <- par()$usr
	ux1 <- uxy[1]; ux2 <- uxy[2]
	uy1 <- uxy[3]; uy2 <- uxy[4]
	px1 <- ux1 + x1 * (ux2 - ux1)
	px2 <- ux1 + x2 * (ux2 - ux1)
	py1 <- uy1 + y1 * (uy2 - uy1)
	py2 <- uy1 + y2 * (uy2 - uy1)
	if(par()$xlog) { px1 <- 10^px1; px2 <- 10^px2 }
	if(par()$ylog) { py1 <- 10^py1; py2 <- 10^py2 }
	arrows(px1, py1, px2, py2, ...)
	invisible(NULL) }
#----------------------------------------addArrows


#addLabel-------------------------------2007-08-22
# Panel label function (Adapted from code by Rob Kronlund)
# Input:
#  x,y - label coordinates in the range (0,1); can step outside
#  txt - desired label at (x,y)
#  ... - arguments used by text, such as "adj", "cex", or "col"
#----------------------------------------RK/JTS/RH
addLabel <- function (x, y, txt, ...) {
	uxy <- par()$usr
	x1 <- uxy[1]; x2 <- uxy[2]
	y1 <- uxy[3]; y2 <- uxy[4]
	x0 <- x1 + x * (x2 - x1)
	y0 <- y1 + y * (y2 - y1)
	if(par()$xlog) x0 <- 10^x0
	if(par()$ylog) y0 <- 10^y0
	text(x0, y0, txt, ...)
	invisible() }
#-----------------------------------------addLabel


#addLegend------------------------------2007-08-22
# Panel key function (Adapted from code by Rob Kronlund)
# Arguments:
#  x,y - label coordinates in the range (0,1); can step outside
#  ... - arguments used by key, such as "lines", "text", or "rectangle"
#----------------------------------------RK/JTS/RH
addLegend <- function (x, y, ...) {
	uxy <- par()$usr
	x1 <- uxy[1]; x2 <- uxy[2]
	y1 <- uxy[3]; y2 <- uxy[4]
	x0 <- x1 + x * (x2 - x1)
	y0 <- y1 + y * (y2 - y1)
	if(par()$xlog) x0 <- 10^x0
	if(par()$ylog) y0 <- 10^y0
	legend(x0, y0, ...)
	invisible(NULL) }
#----------------------------------------addLegend


#drawBars-------------------------------2008-08-16
#  Draw a linear barplot on the current graph
#  x,y   - data coordintates
#  width - bar width, computed if missing
#  base  - y value of the base of each bar
#  ...   - additional parameters for 'lines'
#----------------------------------------------JTS
drawBars <- function (x, y, width, base = 0, ...) {
	nx <- length(x)
	n5 <- 5 * nx
	if ((nx != length(y) || nx == 0)) 
		stop("Inconsistent (x,y)-data.")
	if (missing(width)) {
		width <- ifelse(nx > 1, 0.8 * (x[2] - x[1]), 1)
	}
	if (length(width) == 1) 
		width <- rep(width, nx)
	if (length(base) == 1) 
		base <- rep(base, nx)
	if ((length(width) != nx) || (length(base) != nx)) 
		stop("Inconsistent width or base data.")
	x1 <- numeric(n5)
	y1 <- numeric(n5)
	dx <- width/2
	k <- seq(1:nx)
	x1[5 * k - 4] <- x - dx
	x1[5 * k - 3] <- x - dx
	x1[5 * k - 2] <- x + dx
	x1[5 * k - 1] <- x + dx
	y1[5 * k - 4] <- base
	y1[5 * k - 3] <- y
	y1[5 * k - 2] <- y
	y1[5 * k - 1] <- base
	x1[5 * k] <- NA
	y1[5 * k] <- NA
	xy <- list(x = x1, y = y1)
	lines(xy, ...) }
#-----------------------------------------drawBars


#expandGraph----------------------------2006-08-16
#  Tweaks values to expand margins for multiple graphs
# Arguments:
#  mar - margin paramater
#  mgp - margin points
#  ... - additional par settings
#----------------------------------------------ACB
expandGraph <- function(mar=c(4,3,1.2,0.5), mgp=c(1.6,.5,0),...) {
	par(mar=mar, mgp=mgp, ...)
	invisible() }
#--------------------------------------expandGraph


#lucent---------------------------------2012-07-26
# Create translucent colour
#-----------------------------------------------SM
lucent <- function(col.pal=1,a=1){
	col.rgb<-col2rgb(col.pal)/255
	rgb(t(col.rgb),alpha=a)
}
#-------------------------------------------lucent


#pickCol--------------------------------2006-08-08
# Display interactive colour picking palette
# Arguments:
#  returnValue - if T, user only selects one colour which is returned
#                if F, intermediate GUI is used to display HEX number
#----------------------------------------------ACB
pickCol <- function(returnValue=TRUE) {
	#simply return the first selected value
	if (returnValue)
		return(tclvalue(.Tcl(paste("tk_chooseColor", .Tcl.args(title="Choose a colour")))))
	#otherwise have an intermediate window to display colour codes in
	tt <- tktoplevel()
	tkwm.title(tt,"pickCol()")
	colour <- "#8cda36"
	entryVar<-tclVar(colour)
	entry <- tkentry(tt,textvariable=entryVar, width=8,bg=colour,fg="#000000")
	.changeColour <- function()
	{
		#launch colour picker
		assign("colour",tclvalue(.Tcl(paste("tk_chooseColor",
			.Tcl.args(initialcolor=colour,title="Choose a colour")))),envir=.PBSmodEnv) #.GlobalEnv)
		tmp <- col2rgb(colour)
		#pick white or black foreground colour
		#255*3/2=382.5
		if (sum(tmp)>382 || tmp[2]>180)
			colourFG <- "#000000"
		else
			colourFG <- "#FFFFFF"

		if (nchar(colour)>0) {
			tkconfigure(entry,bg=colour,fg=colourFG)
			tclvalue(entryVar) <- colour
		}
	}
	button <- tkbutton(tt,text="Pick Colour",command=.changeColour)
	tkgrid(entry,button) }
#------------------------------------------pickCol


#plotACF--------------------------------2009-02-24
# Plot auto correlations (support for BRugs)
#-----------------------------------------------RH
plotACF <- function(file,lags=20,clrs=c("blue","red","green","magenta","navy"),...) {
	if (is.vector(file)) file <- matrix(file,ncol=1);
	nc   <- ncol(file); nch <- nc; nr <- nrow(file); lags <- min(nr-1,lags);
	clrs=rep(clrs,nc)[1:nc]
	clim <- qnorm(c(.025,.975))/sqrt(nr);
	acfout <- acf(file,lag.max=lags,plot=FALSE); acfacf <- acfout$acf
	ymin <- min(diag(apply(acfacf,2:3,min)),-.2); ymax <- max(diag(apply(acfacf,2:3,max)));
	xlim <- c(0,lags+.5); ylim <- c(ymin,ymax); ylim[1] <- min(ylim[1],clim[1]); ylim[2] <- max(ylim[2],clim[2]);
	evalCall(plot,argu=list(x=0,y=0,xlim=xlim,ylim=ylim,type="n",tck=.03,xlab="",ylab="",las=1),...,checkdef=TRUE,checkpar=TRUE)
	#plot(0,0,xlim=xlim,ylim=ylim,type="n",tck=.03,xlab="",ylab="",las=1,...)
	if (lags<=30) axis(1,at=1:30,tcl=.25,labels=FALSE);
	abline(h=clim,col="#400080",lty=2);
	for (i in 1:nc) {
		x <- (0:lags)+(i-1)*(.7/nch); y <- acfacf[,i,i];
		evalCall(lines,argu=list(x=x,y=y,type="h",col=clrs[i]),...,checkdef=TRUE,checkpar=TRUE); };
		#lines(x,y,type="h",col=clrs[i],...); };
   abline(h=0,col="grey40",lty=3); box(); };
#------------------------------------------plotACF


#plotAsp--------------------------------2008-08-16
# Plots x and y vectors with plot() but maintaining a fixed aspect
# Arguments:
#  x   - the x coordinates of points in the plot
#  y   - the y coordinates of points in the plot
#  asp - the y/x aspect ratio
#  ... - any arguments to be passed to plot()
#----------------------------------------------ACB
plotAsp <- function(x,y,asp=1,...) {
	dots <- list(...)
	if (is.null(dots$xlim))
		dots$xlim = range(x)
	if (is.null(dots$ylim))
		dots$ylim = range(y)
	xAxisSize <- abs(dots$xlim[1] - dots$xlim[2])
	yAxisSize <- abs(dots$ylim[1] - dots$ylim[2])
	#leave some room for margins
	width <- par("pin")[1]
	height <- par("pin")[2]
	if (xAxisSize > asp*yAxisSize) {
		#x larger than y
		fact <- xAxisSize/(asp*yAxisSize)
		if (width/fact > height) {
			width <- height * fact
		}
		newMaiTop <- (par("fin")[2] - width/fact)
		newMaiSide <- (par("fin")[1] - width)
	}
	else {
		#y larger than x
		fact <- (asp*yAxisSize) / xAxisSize
		if (height/fact > width) {
			height <- width * fact
		}
		#par(pin=c(height/fact, height))
		newMaiTop <- (par("fin")[2] - height)
		newMaiSide <- (par("fin")[1] - height/fact)
	}
	old_mai <- par()$mai
	par(mai=c(par()$mai[1] + (newMaiTop -par()$mai[1]-par()$mai[3])/2,
	          par()$mai[2] + (newMaiSide-par()$mai[2]-par()$mai[4])/2,
	          par()$mai[3] + (newMaiTop -par()$mai[1]-par()$mai[3])/2,
	          par()$mai[4] + (newMaiSide-par()$mai[2]-par()$mai[4])/2))
	plot(x,y,asp=asp,...) 
	par(mai=old_mai) 
	invisible() }
#------------------------------------------plotAsp


#plotBubbles----------------------------2010-08-04
# Function to construct a bubble plot for a matrix z
#  z:     input matrix or data frame
#  xval:  x-values for the columns of z (if length xval != # colums in z, xval is ignored)
#         if xval=TRUE, first row contains x-values for the columns
#  yval:  y-values for the rows of z (if length yval != # rows in z, yval is ignored)
#         if yval=TRUE, first column contains y-values for the rows
#  dnam:  if TRUE, use dimnames as xval and yval
#         (overwrites previously specified values)
#  rpro:  if rpro=TRUE, convert rows to proportions
#  cpro:  if cpro=TRUE, convert columns to proportions
#  rres:  if rres=TRUE, use row residuals (subtract row means)
#  cres:  if cres=TRUE, use column residuals (subtract column means)
#  powr:  power tranform; radii proportional to z^powr
#         powr=0.5 gives bubble areas proportional to z
#  clrs:  colours used for positive, negative, and zero values
#  size:  size (inches) of the largest & smallest bubble
#  lwd:   line width for drawing circles
#  hide0: if TRUE, hide zero-value bubbles
#  frange: fraction by which the range of the axes should be extended
#  prettyaxis: logical: if TRUE, apply the pretty function to both axes
#  ...:   further parameters for the plotting functions 
#-------------------------------------------JTS/RH
plotBubbles <- function(z, xval=FALSE, yval=FALSE, dnam=FALSE, 
    rpro=FALSE, cpro=FALSE, rres=FALSE, cres=FALSE, 
    powr=0.5, size=0.2, lwd=1, clrs=c("black","red","blue"), 
    hide0=FALSE, frange=0.1, prettyaxis=FALSE, ...) {

	if (is.data.frame(z)) {
		use = !sapply(z,is.factor) & sapply(z,is.numeric)
		z=z[,use,drop=FALSE]; if (ncol(z)==0) {showAlert("data frame not useable"); return()}
		#names(z)=gsub("\\.","_",names(z))
		z=as.matrix(z) }
	dz <- dim(z);  ny=ny1=dz[1];  nx=nx1=dz[2]
	if (length(dz)>2) {showAlert("Input matrix must have only 2 dimensions"); return() }
	xval1 <- 1:nx;  yval1 <- 1:ny
	# If first row contains x-values for columns
	if (mode(xval) == "logical") {
		if (xval[1]) {
			xval1 <- z[1,]; ny1 <- ny - 1; } }
	# If first column contains y-values for rows 
	if (mode(yval) == "logical") {
		if (yval[1]) {
			yval1 <- z[,1]; nx1 <- nx - 1; } }
	xind <- (nx - nx1 + 1):nx
	x2=xlabel=xval1[xind]
	yind <- (ny - ny1 + 1):ny
	y2=ylabel=yval1[yind]
	if ((mode(xval) != "logical") & (length(xval) == nx1)) {
		if (mode(xval)=="numeric") x2=xval
		xlabel=xval }
	if ((mode(yval) != "logical") & (length(yval) == ny1)) {
		if (mode(yval)=="numeric") y2=yval
		ylabel=yval }

	zz <- array(z[yind,xind],dim=c(length(yind),length(xind)),dimnames=dimnames(z[yind,xind,drop=FALSE]))
	dots=list(...)
	xlab=dots$xlab; if (is.null(xlab)) xlab=""  # x-axis label
	ylab=dots$ylab; if (is.null(ylab)) ylab=""  # y-axis label

	# dimnames are to be used to over-ride xval and yval
	if (dnam & !is.null(dimnames(zz))) { 
		warn=options()$warn; options(warn=-1)
		if (!is.null(dimnames(zz)[[2]])) {
			xpos = try(as.numeric(dimnames(zz)[[2]]),silent=TRUE)
			if (all(is.na(xpos))) xlabel=dimnames(zz)[[2]]
			else if (!any(is.na(xpos)) && all(diff(xpos)>0 | all(diff(xpos)<0))) {
				xlabel=as.character(xpos); x2=xpos } # strictly increasing / decreasing
		}
		if (!is.null(dimnames(zz)[[1]])) { 
			ypos = try(as.numeric(dimnames(zz)[[1]]),silent=TRUE)
			if (all(is.na(ypos))) ylabel=dimnames(zz)[[2]]
			else if (!any(is.na(ypos)) && all(diff(ypos)>0 | all(diff(ypos)<0))) {
				ylabel=as.character(ypos); y2=ypos } # strictly increasing / decreasing
		}
		options(warn=warn)
	}
#browser()
	xx <- rep(x2, each = length(y2))
	yy <- rep(y2, length(x2))
	minz <- min(zz,na.rm=TRUE);  maxz <- max(zz,na.rm=TRUE);
	if (rpro | cpro) {
		if (minz < 0) {
			zz <- zz - minz
			minz <- 0
			maxz <- max(zz,na.rm=TRUE) } }
	if (rpro) {
		zs <- apply(zz, 1, sum, na.rm=TRUE)
		zz <- sweep(zz, 1, zs, "/") }
	if (cpro) {
		zs <- apply(zz, 2, sum, na.rm=TRUE)
		zz <- sweep(zz, 2, zs, "/") }
	if (rres) {
		zm <- apply(zz, 1, mean, na.rm=TRUE)
		zz <- sweep(zz, 1, zm, "-") }
	if (cres) {
		zm <- apply(zz, 2, mean, na.rm=TRUE)
		zz <- sweep(zz, 2, zm, "-") }
	zNA <- is.na(zz) | is.nan(zz) | is.infinite(zz); zz[zNA] <- 0;
	z0 <- sign(zz) * abs(zz)^abs(powr)
	z1 <- z3 <- z0;  z1[z0 <= 0] <- NA; z3[z0<0 | z0>0] <- NA;
	z2 <- -z0; z2[z0 >= 0] <- NA;
	za <- max(z0,na.rm=TRUE);  zb <- min(z0,na.rm=TRUE)
	zM <- max(abs(z0))
	sz1 <- max(za * size/zM, 0.001)
	sz2 <- max(-zb * size/zM, 0.001)
	#plot(0,0,xlim=extendrange(x2),ylim=extendrange(y2),type="n",xaxt="n",...)
	#axis(1,at=x2,labels=xlabel,...)
	#symbols(xx,yy,circles=as.vector(abs(z0)),inches=size,fg=0,...)
	evalCall(plot,argu=list(x=0,y=0,xlim=extendrange(x2,f=frange),
		ylim=extendrange(y2,f=frange),type="n",axes=FALSE,xlab=xlab,ylab=ylab),...,checkdef=TRUE,checkpar=TRUE)
	if (prettyaxis) {
		if (length(min(x2):max(x2))<=5) xshow = is.element(x2,x2)
		else                            xshow = is.element(x2,pretty(x2,n=10))
		yshow = is.element(y2,pretty(y2,n=10))
	}
	else {
		xshow = rep(TRUE,length(x2)); yshow = rep(TRUE,length(y2)) }
	if (!all(xshow))
		axis(1,at=x2[!xshow],labels=FALSE,tcl=ifelse(is.null(dots$tcl),par()$tcl,dots$tcl)/3)
	if (!all(yshow))
		axis(2,at=y2[!yshow],labels=FALSE,tcl=ifelse(is.null(dots$tcl),par()$tcl,dots$tcl)/3)
	evalCall(axis,argu=list(side=1,at=x2[xshow],labels=xlabel[xshow]),...,checkpar=TRUE)
	evalCall(axis,argu=list(side=2,at=y2[yshow],labels=ylabel[yshow]),...,checkpar=TRUE)
	if (!hide0 && !all(is.na(z3))) {
		evalCall(symbols,argu=list(x=xx,y=yy,circles=as.vector(z3),inches=0.001,fg=clrs[3],lwd=lwd,add=TRUE),...,checkpar=TRUE) }
		#symbols(xx, yy, circles = as.vector(z3), inches = 0.001, fg = clrs[3], lwd = lwd, add = TRUE, ...) }
	if (!all(is.na(z2))) {
		evalCall(symbols,argu=list(x=xx,y=yy,circles=as.vector(z2),inches=sz2,fg=clrs[2],lwd=lwd,add=TRUE),...,checkpar=TRUE) }
		#symbols(xx, yy, circles = as.vector(z2), inches = sz2, fg = clrs[2], lwd = lwd, add = TRUE, ...) }
	if (!all(is.na(z1))) {
		evalCall(symbols,argu=list(x=xx,y=yy,circles=as.vector(z1),inches=sz1,fg=clrs[1],lwd=lwd,add=TRUE),...,checkpar=TRUE) }
		#symbols(xx, yy, circles = as.vector(z1), inches = sz1, fg = clrs[1], lwd = lwd, add = TRUE, ...) }
	box(); invisible(z0) }
#--------------------------------------plotBubbles


#plotCsum-------------------------------2006-07-26
# Plots cumulative frequecy of data
# Arguments:
#  x    - vector of values
#  add  - if TRUE, add cumul. frequency curve to current plot
#  ylim - limits for y-axis
#  xlab - label for x-axis
#  ylab - label for y-axis
#-------------------------------------------RH/ACB
plotCsum<-function(x,add=FALSE,ylim=c(0,1),xlab="Measure",ylab="Cumulative Proportion",...) {
	x <- sort(x); n <- length(x); y <- (1:n)/n
	z <- y >= ylim[1] & y <= ylim[2]
	mdx <- median(x, na.rm = TRUE)
	mnx <- mean(x, na.rm = TRUE); mny <- approx(x,y,xout=mnx)$y
	if (!add) {
		resetGraph();
		evalCall(plot,argu=list(x=x[z],y=y[z],type="n",xlab="",ylab="",las=1,mgp=c(0,.6,0)),...,checkdef=TRUE,checkpar=TRUE)
		#plot(x[z], y[z], type = "n", xlab = "", ylab = "", las=1, mgp=c(0,.6,0), ...)
	}
	evalCall(lines,argu=list(x=x[z],y=y[z],col="blue"),...,checkdef=TRUE,checkpar=TRUE)
	#lines(x[z], y[z], col = "blue")
	abline(h = c(0.5,mny), lty = 3, col=1:2)
	abline(v = c(mdx,mnx), lty = 2, col=1:2)
	addLabel(0.95,0.1,paste("Median = (",paste(signif(c(mdx,.5),3),collapse=", "),")"),cex=1,adj=1)
	addLabel(0.95,0.05,paste("Mean = (",paste(signif(c(mnx,mny),3),collapse=", "),")"),cex=1,adj=1,col=2)
	mtext(xlab, side = 1, line = 2.75, cex = 1.5);  mtext(ylab, side = 2, line = 2.5, cex = 1.5)
	invisible(data.frame(x=x,y=y)) }
#-----------------------------------------plotCsum


#plotDens-------------------------------2009-02-24
# Plot density curves (support for BRugs)
#-----------------------------------------------RH
plotDens <- function(file,clrs=c("blue","red","green","magenta","navy"),...) {
	if (is.vector(file)) file <- matrix(file,ncol=1);
	nc <- ncol(file)
	clrs=rep(clrs,nc)[1:nc]
	dd <- density(unlist(file[,1:nc]),adjust=1.25)
	xlim <- range(dd$x,na.rm=TRUE); ylim <- range(dd$y,na.rm=TRUE)
	for (i in 1:nc) {
		d <- density(file[,i], adjust=1.25);
		xlim[1] <- min(xlim[1],min(d$x)); xlim[2] <- max(xlim[2],max(d$x));
		ylim[1] <- min(ylim[1],min(d$y)); ylim[2] <- max(ylim[2],max(d$y)); };
	evalCall(plot,argu=list(x=0,y=0,xlim=xlim,ylim=ylim,type="n",tck=.03,xlab="",ylab="",las=1),...,checkdef=TRUE,checkpar=TRUE)
	evalCall(lines,argu=list(x=dd$x,y=dd$y,col="grey",lwd=2),...,checkdef=TRUE,checkpar=TRUE)
	#plot(0,0,xlim=xlim,ylim=ylim,type="n",tck=.03,xlab="",ylab="",las=1,...)
	#lines(dd$x,dd$y,col="grey",lwd=2,...)
	for (i in 1:nc) {
		y <- file[,i]; d <- density(y, adjust=1.25)
		evalCall(lines,argu=list(x=d$x,y=d$y,col=clrs[i]),...,checkdef=TRUE,checkpar=TRUE) }; 
		#lines(d$x,d$y,col=clrs[i],...) }; 
	invisible() }
#-----------------------------------------plotDens


#plotFriedEggs--------------------------2015-01-20
#  Pairs plot featuring fried eggs and beer.
#  Original code by Steve Martell (UBC).
#--------------------------------------------SM/RH
plotFriedEggs <- function(A, eggs=TRUE, rings=TRUE,
		levs=c(0.01,0.1,0.5,0.75,0.95), pepper=200, replace=FALSE,
		jitt=c(1,1), bw=25, histclr=NULL)
{
	if (!requireNamespace("KernSmooth", quietly=TRUE))
		stop("Package `KernSmooth' is needed for this function'")
	expandGraph(las=1,mgp=c(0,.75,0))

	panel.cor <- function(x, y, digits=2, prefix="", cex.cor) {
		usr <- par("usr"); on.exit(par(usr))
		par(usr = c(0, 1, 0, 1))
		r <- (cor(x, y))
		txt <- format(c(r, 0.123456789), digits=digits)[1]
		txt <- paste(prefix, txt, sep="")
		if(missing(cex.cor)) cex <- 0.8/strwidth(txt)
		#text(0.5, 0.5, txt, cex = cex * r)
		beer(round(r,2)) }

	panel.hist <- function(x) {
		usr <- par("usr"); on.exit(par(usr))
		par(usr = c(usr[1:2], 0, 1.25))
		h <- hist(x,breaks=20,plot=FALSE)
		breaks <- h$breaks; nB <- length(breaks)
		y <- h$counts; y <- y/max(y)
		if (!is.null(histclr)) {clrs=histclr; bord=1}
		else if (eggs) {clrs=c("moccasin","burlywood"); bord="saddlebrown" }
		else if (rings) {clrs=c("lightsteelblue1","steelblue"); bord="darkblue" }
		else {clrs=c("grey85","grey40"); bord="black" }
		rect(breaks[-nB],0,breaks[-1],y,col=clrs,border=bord);  box() }

	fried.eggs <- function(x, y) {
		bwx=(max(x)-min(x))/bw; bwy=(max(y)-min(y))/bw
		est <- KernSmooth::bkde2D(cbind(x,y),bandwidth=c(bwx,bwy),gridsize=c(51, 51))
		est$fhat=est$fhat/max(est$fhat)
		levs=sort(levs); maxct=max(levs)
		nlev=length(levs); is.white=rev(is.element(nlev:1,seq(nlev,1,-2))); is.yolk=!is.white
		thelines=contourLines(est$x1,est$x2,est$fhat,levels=levs)

		crap=colorRamp(c("sandybrown","gold","white"),space="Lab")
		yolks=apply(crap(seq(.3,.8,len=nlev)),1,function(x){rgb(x[1],x[2],x[3],maxColorValue=255)})
		crap=colorRamp(c("snow","whitesmoke"),space="Lab")
		whites=apply(crap(seq(0,1,len=nlev)),1,function(x){rgb(x[1],x[2],x[3],maxColorValue=255)})

		for (i in 1:(nlev-1)) {
			ii=nlev-i+1
			if (is.yolk[i]) clr=yolks else clr=whites
			polygon(thelines[[i]]$x,thelines[[i]]$y,col=clr[ii],border="grey",lwd=1) }
		if(pepper>0) pepper.mill(x,y,scale=0:7/10)
		polygon(thelines[[nlev]]$x,thelines[[nlev]]$y,col="white",border="grey",lwd=1)
		box() }

	smoke.rings <- function(x, y) {
		bwx=(max(x)-min(x))/bw; bwy=(max(y)-min(y))/bw
		est <- KernSmooth::bkde2D(cbind(x,y),bandwidth=c(bwx,bwy),gridsize=c(51, 51))
		est$fhat=est$fhat/max(est$fhat)
		levs=sort(levs); maxct=max(levs)
		nlev=length(levs)
		points(x,y,pch=20,col="grey",cex=0.2)
		crap=colorRamp(c("white","cornflowerblue","black"),space="Lab") # smoke ring ramp function bounded by white and black
		clrs=apply(crap(seq(.2,.7,len=nlev)),1,function(x){rgb(x[1],x[2],x[3],maxColorValue=255)})
		if (pepper>0) pepper.mill(x,y,scale=0)
		contour(est$x1,est$x2,est$fhat,drawlabels=FALSE,add=TRUE,levels=levs,lwd=2,col=clrs)
		box() }

	pepper.mill = function (x,y,scale=0:9/10) { # grind some pepper
		N=length(x)
		xi=sample(1:N,pepper,replace=ifelse(pepper>N,TRUE,replace))
		points(jitter(x[xi],factor=jitt[1]),jitter(y[xi],factor=ifelse(is.na(jitt[2]),jitt[1],jitt[2])),
			pch=20,cex=0.5,col=grey(scale)) }

	beer <- function(v, col.glass=c("#C6D5D8","#6D939E")) {
		usr=par("usr")
		ymin=usr[3]; ymax=usr[4]; yrng=ymax-ymin
		ymin1=usr[3]+.1*yrng
		ymax1=usr[4]-.2*yrng
		yrng1=ymax1-ymin1
		ymax2=ymin1+abs(v)*yrng1
		xmid=(usr[2]+usr[1])/2
		ymid=(ymax1+ymin1)/2
		xrng=(usr[2]-usr[1])
		Lbot=xmid-.15*xrng; Rbot=xmid+.15*xrng; Ltop=xmid-.25*xrng; Rtop=xmid+.25*xrng             # sides of glass
		Lbeer=xmid-(.15+abs(v)*.1)*xrng;  Rbeer=xmid+(.15+abs(v)*.1)*xrng                          # top of beer

		curved=function(x,yval,scale=0.2) {
			xmin=min(x,na.rm=TRUE); xmax=max(x,na.rm=TRUE)
			s=(2/pi)*asin(sqrt((x-xmin)/(xmax-xmin))) # scale between 0 and 1
			smid=mean(s,na.rm=TRUE)
			ycur=scale/cos(s-smid); ycur=yval+ycur-max(ycur,na.rm=TRUE); return(ycur) }

		xbot=seq(Lbot,Rbot,len=100); ybot=curved(xbot,ymin1) # bottom of glass
		xtop=seq(Ltop,Rtop,len=100); ytop=curved(xtop,ymax1) # top of glass
		xglass=c(xbot,rev(xtop)); yglass=c(ybot,rev(ytop))   # curved glass poly

		#xbeer=c(Lbot,Lbeer,Rbeer,Rbot,Lbot); ybeer=c(ymin1,ymax2,ymax2,ymin1,ymin1) # rectangular beer poly
		xsip=seq(Lbeer,Rbeer,len=100); 
		if (v==0) ysip=ybot else ysip=curved(xsip,ymax2)  # top of beer
		xbeer=c(xbot,rev(xsip)); ybeer=c(ybot,rev(ysip))        # curved beer poly
		bubblex=runif(round(500*abs(v),0),xmid-(.15+abs(v*.95)*.1)*xrng,xmid+(.15+abs(v*.95)*.1)*xrng)
		if (v==0) yfroth=ybot else yfroth=curved(bubblex,ymax2)
		#bubbley=runif(round(500*abs(v),0),ymax2-.02*yrng1,ymax2+.02*yrng1)
		bubbley=runif(round(500*abs(v),0),yfroth-.02*yrng1,yfroth+.02*yrng1)

		polygon(xglass,yglass,,col="aliceblue",border=FALSE)              # glass surface
		lines(xtop,ymax1+abs(ytop-ymax1),lwd=2,col=col.glass[1])          # top back rim
		if (v!=0) {
			polygon(xbeer,ybeer,col=ifelse(v<0,"orange","gold"),border="burlywood")                # beer
			points(bubblex,bubbley,pch=21,col=ifelse(v<0,"orange","gold"),bg="white",cex=seq(0.1,1,length=10)) }
		lines(xbot,ybot,lwd=3,col=col.glass[1])                           # bottom edge
		lines(xbot,ybot-.005,lwd=1,col=col.glass[2])                      # bottom edge shadow
		lines(c(Lbot,Ltop),c(ymin1,ymax1),lwd=4,col=col.glass[1])         # left edge
		lines(c(Lbot-.01,Ltop-.01),c(ymin1,ymax1),lwd=1,col=col.glass[2]) # left edge shadow
		lines(c(Rbot,Rtop),c(ymin1,ymax1),lwd=4,col=col.glass[1])         # right edge
		lines(c(Rbot+.01,Rtop+.01),c(ymin1,ymax1),lwd=1,col=col.glass[2]) # right edge shadow
		lines(xtop,ytop,lwd=3,col=col.glass[1])                           # top front rim
		text(xmid,ymid,labels=c(paste(ifelse(v<0,"-","+"),abs(v))),cex=abs(v*100)^.1)
	}
	if (eggs) lower=fried.eggs else if (rings) lower=smoke.rings else lower=pepper.mill
	pairs(A,lower.panel=lower,diag.panel=panel.hist,upper.panel=panel.cor,gap=0,label.pos=0.92)
}
#------------------------------------plotFriedEggs


#plotSidebars---------------------------2012-07-26
# Plot (x,y) table as horizontal sidebars.
# Additional arguments:
#   lbl    - labels for the x- and y-axis
#   margin - function to report margin summaries
#--------------------------------------------SM/RH
plotSidebars <- function(z, scale=1, col=lucent("blue",0.25), ...){
	ciao = function(){ par(oldpar) }
	oldpar = par(no.readonly=TRUE)
	on.exit(ciao())
	dots = list(...)
	plot.new()
	# set up graphical parameters
	#par(mar=c(1, 1, 1, 1), oma=c(1, 1, 1, 1)*4, xaxs='r')
	evalCall(par,argu=list(no.readonly=FALSE),...,checkpar=TRUE)
	dz  <- dim(z)
	dn  <- dimnames(z)
	if(!is.null(dots$lbl)) lbl = dots$lbl
	else {
		lbl <- rev(names(dn))
		if(is.null(lbl)) lbl=c("Year","Age") }
	
	ny = nr = dz[1]
	nx = nc = dz[2]
	y  = as.numeric(unlist(dn[1]))
	x  = as.numeric(unlist(dn[2]))
	yl = extendrange(y, f=0.05)
	xl = extendrange(x, f=0.05*sqrt(scale))
	
	par(usr=c(xl, yl))
	box()
	mtext(lbl, side=1:2, outer=FALSE, line=2.5, cex=1.5, las=0)
	axis(1); axis(2)
	scale = scale/diff(range(z,na.rm=TRUE))
	for( i in 1:nx ) {
		dy = 0.5*diff(y[1:2])
		xx = as.vector(rbind(y-dy, y+dy))
		xx = c(min(xx), xx, max(xx))
		yy = as.vector(rbind(z[,i], z[,i]))*scale
		yy = c(0, yy, 0)
		#polygon(x[i]-yy, xx,col=col)
		evalCall(polygon,argu=list(x=x[i]-yy, y=xx, density=NULL,angle=45,border=NULL,col=col,lty=par("lty"),fillOddEven=FALSE),
			..., checkpar=TRUE)
	}
	if (!is.null(dots$margin) && is.function(dots$margin)) {
		xmarg = apply(z,2,dots$margin)
		ymarg = apply(z,1,dots$margin)
		usr=par()$usr
		text(usr[2]-0.02*diff(usr[1:2]), y, ymarg,adj=1,cex=0.7,font=2)
		text(x, usr[3]+0.02*diff(usr[3:4]), xmarg,adj=1,cex=0.7,font=2)
	}
}
#-------------------------------------plotSidebars


#plotTrace------------------------------2009-02-24
# Plot trace lines (support for BRugs)
#-----------------------------------------------RH
plotTrace <- function(file,clrs=c("blue","red","green","magenta","navy"),...) {
	if (is.vector(file)) file <- data.frame(x=1:length(file),y=file);
	nc <- ncol(file)
	clrs=rep(clrs,nc)[1:nc]
	x  <- file[,1]; xlim <- range(x); ylim <- range(file[,2:nc])
	evalCall(plot,argu=list(x=0,y=0,xlim=xlim,ylim=ylim,type="n",tck=.03,xlab="",ylab="",las=1),...,checkdef=TRUE,checkpar=TRUE)
	#plot(0,0,xlim=xlim,ylim=ylim,type="n",tck=.03,xlab="",ylab="",las=1,...)
	for (i in 2:nc) {
		y <- file[,i]
		evalCall(lines,argu=list(x=x,y=y,col=clrs[i-1]),...,checkdef=TRUE,checkpar=TRUE) }
		#lines(x,y,col=clrs[i-1],...) }
	invisible() }
#----------------------------------------plotTrace


#resetGraph-----------------------------2012-12-04
# Resets par() values to R default
#-------------------------------------------ACB/RH
resetGraph=function(reset.mf=TRUE)
{
	#ensure init has been called (to pass R check)
	#PBSmodelling:::.initPBSoptions()
	.initPBSoptions()
	tget(.PBSmod)

	#cache value on first run
	if (is.null(.PBSmod[[".options"]][["par.default"]])) {
		dev.new()
		p <- graphics::par(no.readonly = TRUE)
		dev.off()
		#eval(parse(text=".PBSmod$.options$par.default <<- p"))
		.PBSmod$.options$par.default <- p
		tput(.PBSmod)
	}
	if (dev.cur()==1) frame()
	else {
		p = .PBSmod$.options$par.default
		if (reset.mf) {
			par(p) ; frame()
		} else {
			keep=c("mfrow","mfcol","mfg") #keep these settings
			#keep=c("mfrow","mfg") #keep only mfrow
			fixed=par(keep)
			reset=p[setdiff(names(p),keep)]
			par(reset); par(fixed); par(new=FALSE)
		}
	}
	invisible(par()) }
#--------------------------------------resetGrapgh


#testAlpha------------------------------2011-09-08
# Display various alpha transparencies
#-----------------------------------------------RH
testAlpha <- function (alpha=seq(0,1,len=25), fg="blue", bg="yellow",
      border="black", grid=FALSE, ...) {
	par0 <- par(no.readonly = TRUE)
	on.exit(par(par0))
	N=length(alpha); fg=rep(fg,N)[1:N]; bg=rep(bg,N)[1:N]; border=rep(border,N)[1:N]
	rc=.findSquare(N); m=rc[1]; n=rc[2]
	rgbfg=col2rgb(fg)/255; rgbbg=col2rgb(bg)/255; rgbbo=col2rgb(border)/255
	rgba=rbind(rgbfg,alpha,rgbbg,rgbbo)
	resetGraph()
	expandGraph(mfrow=rc,mar=c(0,0,0,0),oma=c(3,3,.5,.5))
	apply(rgba,2,function(x){
		evalCall(plot,argu=list(x=0,y=0,xlim=c(0,1),ylim=c(0,1),type="n",xlab="",ylab="",axes=FALSE),...,
			checkdef=TRUE,checkpar=TRUE)
		evalCall(polygon,argu=list(x=c(.1,.9,.9,.1),y=c(.1,.1,.9,.9), 
			col=rgb(x[5],x[6],x[7]), border=FALSE), ..., checkpar=TRUE)
		z1=rep(.1,9); z2=seq(.1,.9,.1); z3=rep(.9,9)
		if (grid) {
			segments(x0=z1,y0=z2,x1=z3,y1=z2,col="grey"); segments(x0=z2,y0=z1,x1=z2,y1=z3,col="grey") }
		evalCall(polygon,argu=list(x=c(.25,.75,.75,.25), y=c(.25,.25,.75,.75),
			col=rgb(x[1],x[2],x[3],x[4]), border=rgb(x[8],x[9],x[10],x[4])), ..., checkpar=TRUE)
		evalCall(addLabel,argu=list(x=.5,y=.05,txt=eval(parse(text=paste("expression(alpha==",round(x[4],3),")"))),
			cex=1.5), ..., checkpar=TRUE)
	})
	invisible(rgba) }
#----------------------------------------testAlpha


#testCol--------------------------------2011-09-08
# Display test colours as circular patches.
#  cnam - colour names to search for
#-----------------------------------------------RH
testCol <- function(cnam=colors()[sample(length(colors()),15)]) {
	#get similar colours
	getCol <- function(x) {
		palette <- colors()
		n <- length(palette)
		z <- NULL
		for (i in x) {
			a <- regexpr(i,palette)
			b <- (1:n)[a>0]
			z <- union(z,b)
		}
		lovely <- palette[z]
		return(lovely)
	}
	par0 <- par(no.readonly = TRUE)
	on.exit(par(par0))
	clrs <- getCol(grep("^[^#0-9]", cnam,value=TRUE))                       # identifies all colours that contain string patterns
	clrs <- c(clrs, 
	          grep("^#[0-9a-f]{6,8}$", cnam, value=TRUE, ignore.case=TRUE), # identifies hexadecimal colours, incl. transparents
	          grep("^[0-9]+$", cnam, value=TRUE, ignore.case=TRUE)          # identifies colours specified by an integer
	          )
#browser()
	# fiddle for mfrow
	N <- length(clrs); din <- par()$din; x <- din[1]; y <- din[2]
	cell <- sqrt(prod(din)/N)
	cols <- ceiling(x/cell); rows <- ceiling(y/cell)
	if (N <= rows*cols-cols) rows <- rows-1

	xlim <- c(1, cols) + c(-.25,.25); ylim <- c(-rows,-1) + c(-.25,.25)
	resetGraph()
	par(mfrow=c(1,1),mai=c(.05,.05,.05,.05))
	plot(0,0,xlim=xlim,ylim=ylim,type="n",axes=FALSE,xlab="",ylab="")
	k <- 0
	for (i in 1:rows) {
		for (j in 1:cols) {
			k <- k+1
			points(j,-i, col=clrs[k], pch=16,cex=5)
			text(j,-i-.04*diff(ylim),clrs[k],cex=.6) } }
	invisible(clrs) }
#------------------------------------------testCol


#testLty--------------------------------2011-09-08
# Display line types available
# Arguments:
#  newframe - if T, clear graphics frame, if F, overlay
#-----------------------------------------------RH
testLty <- function (newframe=TRUE, n=1:18, ...) {
	if (newframe) frame()
	par0 <- par(no.readonly = TRUE)
	on.exit(par(par0))
	n = sort(unique(n))
	xusr = c(n[1], rev(n)[1])
	if (length(n)==1) xusr = xusr + c(-0.5,0.5)
	par(usr = c(xusr, 0, 1) )
	for (i in n) lines(c(i, i), c(0, 1), lty = i, ...)
	mtext(as.character(n), side = 1, line = 1, at = n)
	mtext("LINE TYPES (lty)", side = 3, line = 2)
	invisible(NULL) }
#------------------------------------------testLty


#testLwd--------------------------------2011-09-08
# Display line widths
# Arguments:
#  lwd      - line widths to test
#  col      - colours to use
#  newframe - if T, use a new graphics frame, if F, overlay
#-----------------------------------------------RH
testLwd <- function (lwd=1:20, col=c("black","blue"), newframe=TRUE) {
	if (newframe) { resetGraph(); frame(); }
	par0 <- par(no.readonly = TRUE)
	on.exit(par(par0))
	xlim <- range(lwd);
	nl <- length(lwd); nc <- length(col); col <- rep(col,ceiling(nl/nc));
	par(usr = c(c(xlim[1]-1, xlim[2]+1), c(0, 1)))
	for (i in lwd) lines(c(i, i), c(0, 1), lty = 1, lwd = i, col=col[i-xlim[1]+1])
	mtext(as.character(lwd), side = 1, line = 1, at = lwd)
	mtext(paste("LINE WIDTHS (",xlim[1],"-",xlim[2],")"), side=3, line=2)
	invisible(NULL) }
#------------------------------------------testLwd


#testPch--------------------------------2011-09-09
# Display plotting symbols or octal strings
# Arguments:
#  pch      - symbols/octals to test
#  ncol     - number of columns to use
#  grid     - logical: if T, display in a grid
#  newframe - if T, use a new graphics frame, if F, overlay
#  octal    - logical: if T, use octal (backslash) strings
#-----------------------------------------------RH
testPch <- function (pch=1:100, ncol=10, grid=TRUE, newframe=TRUE, octal=FALSE, ...) {
	if (!is.element(ncol,c(2,5,10))) stop("Set ncol to 2 or 5 or 10")
	if (!all(diff(pch)==1)) stop("pch vector must be a continuous increasing integer series")
	if (!octal && (all(pch>255) | any(pch<0))) stop("pch must be in the range 0 - 255")
	if (octal && (all(pch<41) | all(pch>377))) stop("pch must be in the range 41 - 377")
	if (newframe) {
		resetGraph(); frame(); }
	par0 <- par(no.readonly = TRUE)
	on.exit(par(par0))
	npch =length(pch)
	xlim = c(.5,ncol+.5)
	rlim = floor((pch[c(1,npch)]-1)/ncol); yval = rlim[1]:rlim[2]
	ylim = rev(rlim); ylim = ylim + c(.5,-.5)
	pchy = pch[is.element(pch,seq(0,1000,ncol))]
	if(length(pchy)<length(yval)) {
		pchy = c(pchy,floor((pchy[length(pchy)]+ncol)/ncol)*ncol)
	}
	ylab = pchy - ncol
	par(usr=c(xlim,ylim))
	if (grid) {
		abline(v=seq(.5,ncol+.5,1),h=seq(rlim[1]-.5,rlim[2]+.5,1),col="gainsboro")
	}
	for (i in pch) {
		y = floor((i - 1)/ncol)
		x = i - ncol * y
		if (octal) {
			#if (i<41 | i>377 | is.element(i,seq(9,379,10)) | is.element(i,c(90:99,190:199,290:299))) next
			if (i<41 | i>377 | grepl("[89]",i) | is.element(i,c(177,201,220,235))) next   # Duncan Murdoch suggested grepl("[89]",i)
			cset = eval(parse(text=paste("\"\\", i, "\"", sep = "")))
			text(x,y, cset, cex=1.5, ...) 
		}
		else {
			if (i>255 | is.element(i,26:31)) next
			points(x, y, pch = i, cex=1.5, ...)
		}
	}
	mtext(as.character(1:ncol), side=1, line=.5, at=(1:ncol), cex=1.3, col="blue")
	mtext(as.character(1:ncol), side=3, line=.4, at=(1:ncol), cex=1.3, col="blue")
	mtext(ylab, side=2, line=1, at=yval, cex=1.3, col="red",las=1)
	mtext(paste(ifelse(octal,"OCTAL STRINGS","PCH CHARACTERS"),"(",pch[1],"-",pch[npch],")"), side=3, line=2.2, cex=1.2)
	invisible(yval) }
#------------------------------------------testPch


#===== THE END ===================================

