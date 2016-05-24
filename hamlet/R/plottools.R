#
# Various plotting functions, especially for mixed and/or matched data
#

# Scatterplot function for mixed type (numerical and categorical) variables with marginal distributions in each scatterplot
mixplot <- function(
	# Data frame or a matrix
	x,
	# Main title
	main = NA,
	# If a matching matrix of size nrow(x) x ncol(x) or a matching vector of length(x) is provided,
	# these are annotated in the dimensions by connecting the observations with lines
	# For a matching matrix different values are given as parameter 'par'
	match,
	func = function(x, y, par) { segments(x0=x[1], y0=x[2], x1=y[1], y1=y[2], col=par) },
	# If legend should be constructed and plotted
	legend = T,
	# Colors for the categories or single observations if categories are not present
	col = palette(),
	# Should lines be drawn to represent one of the variables if the other one is missing in a 2-dim scatterplot
	na.lines = T,
	# Plot origin {0,0} using vertical and horizontal ablines
	origin = F,
	# Should marginal distributions be plotted, value of 0/FALSE/"no"/"none", 1/TRUE/"rug", 2/"hist"
	marginal = F,
	# Layout heights
	lhei,
	# Layout widths
	lwid,
	# Level of verbosity: -1<= (no verbosity), 0/FALSE (warnings) or >=1/TRUE (additional information)
	verb = 0,
	# Additional parameters
	...
){
	# Old par-settings (taken from ?layout)
	def.par <- par(no.readonly = TRUE) # save default
	
	# Make sure that some of the parameters are of correct class or cast to one
	verb <- as.numeric(verb)
	marginal <- as.character(marginal)
	if(is.vector(x)) x <- as.matrix(x)

	# Which fields to consider numeric or categorical (treated differently)
	nums <- c("numeric", "integer", "logical")
	ctgr <- c("factor", "character")
	
	# Extract types for the columns, cast x to data.frame if necessary
	if(class(x)=="data.frame"){
		colclass <- lapply(x, FUN=class)
	}else if(class(x)=="matrix"){
		colclass <- class(x[1,1])
	}else{
		x <- as.data.frame(x)
		colclass <- lapply(x, FUN=class)
	}

	# Number of matrix elements in the plotting region ncol_numeric x ncol_numeric
	wnum <- colclass %in% nums
	nnum <- sum(wnum)
	# Number of different categorical variables
	wcat <- colclass %in% ctgr
	ncat <- sum(wcat)

	# If no categorical variables available, introduce artificial category "Observation" which affects all
	if(ncat==0){
		cats <- rep("Observation", times=nrow(x))
		# Old, wrong!
		#wcats <- 1:nrow(x)
		wcats <- rep(1, times=nrow(x))
	}else{
	# Else map observations to categories based on ctgr-type columns
		cats <- apply(x[,wcat,drop=F], MARGIN=1, FUN=function(x) paste(x, collapse=", "))
		wcats <- base::match(cats, unique(cats))
	}
	# Colors for the categories per each observation
	col <- rep(col, length.out=length(wcats))[wcats]

	# If user has not defined alternative lhei (layout heights) and lwid (layout widths) construct them
	if(missing(lhei)){
		lhei <- c()
		if(!is.na(main)) lhei <- c(lhei, 0.1)
		if(!as.character(marginal) %in% c("0", "FALSE", "no", "none")){
			lhei <- c(lhei, rep(c(0.05,0.3), times=ifelse(nnum>2, nnum, 1)))
		}else{
			lhei <- c(lhei, rep(0.3, times=ifelse(nnum>2, nnum, 1)))
		}
		if(legend) lhei <- c(lhei, 0.1)
	}
	if(missing(lwid)){
		if(!as.character(marginal) %in% c("0", "FALSE", "no", "none")){
			lwid <- rep(c(0.3, 0.05), times=ifelse(nnum>2, nnum, 1))
		}else{
			lwid <- rep(0.3, times=ifelse(nnum>2, nnum, 1))
		}
	}

	# Division of device plot region to subregions
	subsize <- ifelse(!as.character(marginal) %in% c("0", "FALSE", "no", "none"), 4, 1)
	if(subsize==4){
		temp <- rbind(rep(letters[1:2], times=ifelse(nnum>2, nnum, 1)), rep(letters[3:4], times=ifelse(nnum>2, nnum, 1)))[rep(1:2, times=ifelse(nnum>2, nnum, 1)),]
		temp[temp=="a"] <- seq(from=1+0, to=4*(ifelse(nnum>2, nnum, 1)*ifelse(nnum>2, nnum, 1))+0, by=4)
		temp[temp=="b"] <- seq(from=1+1, to=4*(ifelse(nnum>2, nnum, 1)*ifelse(nnum>2, nnum, 1))+1, by=4)
		temp[temp=="c"] <- seq(from=1+2, to=4*(ifelse(nnum>2, nnum, 1)*ifelse(nnum>2, nnum, 1))+2, by=4)
		temp[temp=="d"] <- seq(from=1+3, to=4*(ifelse(nnum>2, nnum, 1)*ifelse(nnum>2, nnum, 1))+3, by=4)
		class(temp) <- "numeric"
		lmat <- rbind( 
				rep(1, times=length(lwid)),
				t(temp)+1,
				rep(1+max(temp)+1, times=length(lwid))
			)
	}else{
		lmat <- rbind( 
				rep(1, times=length(lwid)),
				matrix(2:(ifelse(nnum>2, nnum, 1)*ifelse(nnum>2, nnum, 1)+1), ncol=nnum, byrow=T),
				rep(1+ifelse(nnum>2, nnum, 1)*ifelse(nnum>2, nnum, 1)+1, times=length(lwid))
			)
	}
	# Remove header row from layout matrix if it is not desired
	if(is.na(main)){
		lmat <- lmat - 1
		lmat <- lmat[-1,]
	}
	# Remove legend row from layout matrix if it is not desired
	if(!legend){
		lmat <- lmat[-nrow(lmat),]
	}
	if(as.numeric(verb)>=1){
		print("lmat"); print(lmat)
		print("lwid"); print(lwid)
		print("lhei"); print(lhei)
	}
	# Set up the layout matrix to the plot device
	l <- layout(lmat, widths=lwid, heights=lhei)

	# If we should plot the title
	if(!is.na(main)){
		par(mar=c(0,0,0,0)); plot.new(); plot.window(xlim=c(-1,1), ylim=c(-1,1))
		text(0,0, main)
	}

	# If user wants to annotate matches
	if(!missing(match)){
		# Process matching matrix to a matching vector
		if(class(match) %in% c("matrix", "data.frame")){
			match <- match.mat2vec(as.matrix(match))
		}
		# Unique matching elements
		mm <- unique(match)
		# Which element belongs to which unique element
		ms <- lapply(mm, function(z) which(match==z))
		# The 2 rows are pairs, columns indicate each pair-combination of observations to connect
		ms <- lapply(ms, function(z) { if(length(z)>1) { combn(z, 2) } else NA })
		if(verb>=1){
			print("Matching information processed:")
			print(ms)
		}
		plotmatch <- T
	}else{
		plotmatch <- F
	}
	# Sub plotting function for each 1x1 or 2x2 subplot depending on marginal type
	subplox <- function(x, marginal, labels = F, plotmatch = F, ...){
		# Should variable labels be shown (and thus larger margins)
		if(labels){
			mar <- c(4,4,1,1)
		}else{
			mar <- c(2,2,1,1)
		}
		# Top marginal distribution
		if(!as.character(marginal) %in% c("0", "FALSE", "no", "none")){
			martemp <- mar; martemp[c(1,3)] <- 0; par(mar=martemp)			
			# Marginals
			plot.new(); plot.window(xlim=grDevices::extendrange(x[,1]), ylim=c(-1,1))
			# Rug
			if(any(as.character(marginal) %in% c("1", "TRUE", "rug"))){
				for(i in 1:nrow(x)) abline(v=x[i,1], col=col[i], lwd=2)
			}
			# Histogram
			if(any(as.character(marginal) %in% c("2", "hist"))){
				uniqs <- unique(wcats)
				hall <- hist(x[,1], plot=F, breaks=seq(from=min(x[,1], na.rm=T), to=max(x[,1], na.rm=T), length.out=50))
				for(i in 1:length(uniqs)){
					h <- hist(x[wcats==uniqs[i],1], plot=F, breaks=seq(from=min(x[,1], na.rm=T), to=max(x[,1], na.rm=T), length.out=50))
					for(j in 1:length(h$density)){
						d <- h$counts[j]/max(hall$counts)
						rect(xleft = h$breaks[j], xright = h$breaks[j+1], ybottom = -1, ytop = unlist(ifelse(d>0, d, -1)), col=unique(col)[i], border=NA)
					}
				}
			}
		}
		# Actual scatterplot
		par(mar=mar)
		plot(x, xlim=extendrange(x[,1]), ylim=extendrange(x[,2]), ...)
		# Matching combinations, pairwise-applied 'func' to data points belonging to same submatch
		if(plotmatch){
			lapply(1:length(ms), FUN=function(z){
				if(is.matrix(ms[[z]]))
				lapply(1:ncol(ms[[z]]), FUN=function(q){
					func(x=as.numeric(x[ms[[z]][1,q],]), y=as.numeric(x[ms[[z]][2,q],]), par=z)
				})
			})
		}
		# Origin lines if desired
		if(origin){
			abline(h=0, col="grey")
			abline(v=0, col="grey")
		}
		# Top-right empty box and right marginal distribution
		if(!as.character(marginal) %in% c("0", "FALSE", "no", "none")){
			# Empty box top-right
			par(mar=c(0,0,0,0))
			plot.new(); plot.window(xlim=c(-1,1), ylim=c(-1,1))
			martemp <- mar; martemp[c(2,4)] <- 0; par(mar=martemp)			
			# Right marginal
			plot.new(); plot.window(xlim=c(-1,1), ylim=extendrange(x[,2]))
			# Rug
			if(any(as.character(marginal) %in% c("1", "TRUE", "rug"))){
				for(i in 1:nrow(x)){ abline(h=x[i,2], col=col[i], lwd=2)}
			}
			# Histogram
			if(any(as.character(marginal) %in% c("2", "hist"))){
				uniqs <- unique(wcats)
				hall <- hist(x[,2], plot=F, breaks=seq(from=min(x[,2], na.rm=T), to=max(x[,2], na.rm=T), length.out=50))
				for(i in 1:length(uniqs)){
					h <- hist(x[wcats==uniqs[i],2], plot=F, breaks=seq(from=min(x[,2], na.rm=T), to=max(x[,2], na.rm=T), length.out=50))
					for(j in 1:length(h$density)){
						d <- h$counts[j]/max(hall$counts)
						rect(ytop = h$breaks[j], ybottom = h$breaks[j+1], xleft = -1, xright = unlist(ifelse(d>0, d, -1)), col=unique(col)[i], border=NA)
					}
				}
			}
		}
	}
	
	# Plot actual scatterplots
	# If only one numeric column handle this exceptionally
	if(nnum==1){
		subplox(x=x[,which(wnum)[1]], marginal=marginal, labels = T, col=col, plotmatch = plotmatch, ...)
	# Else multiple scatterplots
	}else if(nnum==2){
		subplox(x=x[,which(wnum)[1:2]], marginal=marginal, labels = T, col=col, plotmatch = plotmatch, ...)
	}else{
		for(row in 1:nnum){
			for(column in 1:nnum){
				if(verb>=1) print(paste("row", row, "column", column))
				if(row==column){
					par(mar=c(0,0,0,0))
					# Empty top marginal
					if(!as.character(marginal) %in% c("0", "FALSE", "no", "none")){
						plot.new(); plot.window(xlim=c(-1,1), ylim=c(-1,1))	
					}
					# Text of the variable
					plot.new()
					plot.window(xlim=c(-1,1), ylim=c(-1,1))
					text(0,0,colnames(x)[which(wnum)[row]])
					# Empty top right box and right marginal
					if(!as.character(marginal) %in% c("0", "FALSE", "no", "none")){
						plot.new(); plot.window(xlim=c(-1,1), ylim=c(-1,1))	
						plot.new(); plot.window(xlim=c(-1,1), ylim=c(-1,1))	
					}
				}else{
					subplox(x=x[,rev(which(wnum)[c(row,column)])], marginal=marginal, col=col, plotmatch = plotmatch, ...)
				}
			}
		}
	}
	
	# Possibly build the legend and plot it	
	if(!as.numeric(legend)==0){
		par(mar=c(0,0,0,0))
		plot.new(); plot.window(xlim=c(-1,1), ylim=c(-1,1))	
		if(verb>=1){
			print(col)
			print(wcats)
			print(cats)
		}
		legend("bottom", col=unique(col)[unique(wcats)], legend=unique(cats)[unique(wcats)], horiz=T, bty="n", ...)
	}
	
	# Restore old par settings
	par(def.par)
	
	# Return invisibly some information about the plot
	invisible(list(x=x, lmat=lmat, lwid=lwid, lhei=lhei))
}

# Plot-region based heatmap
# 
# Function tries to compatible with such functions as 'heatmap' in base package or 'heatmap.2' in gplots-package
hmap <- function(
	# Input data matrix to plot
	x,

	# Plotting region settings
	#
	# Whether we want to add to an existing region or create a new one
	add = F,
	# x and y axis limits in the heatmap itself
	# NOTE! This does not include the top and left dendrograms or row and column labels 
	# (see parameters leftlim, toplim, rightlim, bottomlim for these)
	xlim=c(0.2,0.8),
	ylim=c(0.2,0.8),
	
	# Colors to use
	col = heat.colors(10),
	
	# Border and line settings for the heatmap bins
	# Should be matrices of equal dimensions to x
	# NOTE! These are not reordered according to the marginal dendrograms, 
	# but instead affect the absolute rectangle positions at the plot
	#
	# Color for borders in the heatmap bins
	border = matrix(NA, nrow=nrow(x), ncol=ncol(x)),
	# Line type for borders in the heatmap bins
	lty = matrix("solid", nrow=nrow(x), ncol=ncol(x)),
	# Line width for borders in the heatmap bins
	lwd = matrix(1, nrow=nrow(x), ncol=ncol(x)),
	
	# Settings / functions for the marginal hierarchical clusterings
	#
	# Hierarchical clustering function
	hclustfun = hclust, 
	# Distance function to provide to hclustfun
	distfun = dist,
	# Reordering function used to reorder rows and columns while preserving legality of clustering
	reorderfun = function(d, w) reorder(d, w),
	# Function used for plotting name labels
	textfun = function(xseq, yseq, labels, type="row", ...){ if(type=="col") par(srt=90); text(x=xseq, y=yseq, labels=labels, ...); if(type=="col") par(srt=0) },
	# Should symmetricity in rows/columns be enforced
	symm = F,
	# Reordering of rows
	Rowv=NULL,
	# Reordering of columns
	Colv=if(symm) Rowv else NULL,
	# x-axis limits for plotting the dendrogram to the left
	leftlim = c(0, 0.2),
	# y-axis limits for plotting the dendrogram to the top
	toplim = c(0.8, 1),
	# x-axis limits for plotting row names to the right
	rightlim = c(0.8, 1),
	# y-axis limits for plotting column names to the bottom
	bottomlim = c(0, 0.2),
	# Type of clustering visualization, by default "rect"=rectangular, can be also "tri"angular
	type = "rect",
	
	# Should data matrix 'x' be scaled according to 'row's, 'column's or 'none'
	# NOTE! While 'heatmap' by default scales by rows, here assumed user does not want scaling by default
	scale = c("none", "row", "column"),
	# Should missing values be removed
	na.rm = T,
	# Number of discrete bins to divide the data to
	nbins = length(col),
	# Value ranges for the discrete binning
	valseq = seq(from=min(x, na.rm=na.rm), to=max(x, na.rm=na.rm), length.out=nbins),	
	# Should row names be plotted (using text-func),
	# if boolean, then rownames(x) will be plotted. If it is a custom vector of length nrow(x), then these names will be plotted instead
	namerows = T,
	# Should column names be plotted (using text-func)
	# if boolean, then colnames(x) will be plotted. If it is a custom vector of length ncol(x), then these names will be plotted instead
	namecols = T,
	
	# Additional parameters
	...
){
	# Number of rows and columns
	nr = nrow(x)
	nc = ncol(x)
	
	# Recycle border color, line type and line width settings if dimensions are not equal to x
	if(!identical(dim(border),dim(x))){
		print("Border not identical")
		border = matrix(border, nrow=nrow(x), ncol=ncol(x))
	}
	if(!identical(dim(lty),dim(x))){
		print("lty not identical")
		lty = matrix(lty, nrow=nrow(x), ncol=ncol(x))
	}
	if(!identical(dim(lwd),dim(x))){
		print("lwd not identical")
		lwd = matrix(lwd, nrow=nrow(x), ncol=ncol(x))
	}
	
	# If requested, create a new plotting region
	if(!add){
		plot.new()
		plot.window(xlim=range(leftlim,xlim,rightlim), ylim=range(ylim,toplim,bottomlim))
	}
	
	# Recursive function for manually plotting dendrograms (hidden)
	plotdend <- function(
		dend, # dendrogram object
		horiz=F, # Horizontal (by default vertical)
		left=0, # How many leaves are to the left and to the right, used to determine center 
		type="rect", # "rect"angular (or else: triangular)
		new=T, # Make new graphics device
		leafpos = seq(from=0, to=attr(dend, "members"), length.out=(attr(dend, "members")+1)), # positions for the leafs
		normleafpos = T, # Should leaf positions be normalized so the ends will match with the centers of respective box elements
		xlim=c(0, attr(dend, "members")),
		ylim=c(0, attr(dend, "height")),
		xscale = abs(xlim[2] - xlim[1])/attr(dend, "members"), # x-axis scaling factors according to chosen coordinates
		yscale = abs(ylim[2] - ylim[1])/attr(dend, "height"), # y-axis scaling factors according to chosen coordinates
		revx = F, # Should x-axis coordinates be reversed
		revy = F, # Should y-axis coordinates be reversed
		...
		){ 
		if(normleafpos) leafpos <- unlist(lapply(1:(length(leafpos)-1), FUN=function(z) mean(c(leafpos[z], leafpos[z+1]))))
		if(new){
			plot.new()
			if(!horiz) plot.window(xlim=xlim, ylim=ylim)
			if(horiz) plot.window(xlim=ylim, ylim=xlim)
		}
		if(!attr(dend, "members")==1 ){
			# Branches into 2:
			leftmost <- leafpos[left+1]
			if(is.null(attr(dend, "midpoint"))){
				mid <- 0
			}else{
				mid <- attr(dend, "midpoint")
			}
			
			# Compute node positions
			leftleftmost <- left
			rightleftmost <- left + attr(dend[[1]], "members") 
			leftmid <- attr(dend[[1]], "midpoint")
			if(is.null(leftmid)) leftmid <- 0
			rightmid <- attr(dend[[2]], "midpoint")
			if(is.null(rightmid)) rightmid <- 0

			y0a <- ylim[1] + (attr(dend, "height")*yscale)
			y0b <- ylim[1] + (attr(dend[[1]], "height")*yscale)
			x0a <- xlim[1] + (leftmost + mid)*xscale
			x0b <- xlim[1] + (leftmost + leftmid)*xscale
			# For right branch
			x1a <- xlim[1] + (leafpos[rightleftmost+1] + rightmid)*xscale
			y1a <- ylim[1] + (attr(dend[[2]], "height")*yscale)
			if(revx){
				x0a <- xlim[2] - x0a
				x0b <- xlim[2] - x0b
				x1a <- xlim[2] - x1a
			}
			if(revy){
				y0a <- ylim[2] + leftlim[1] - y0a
				y0b <- ylim[2] + leftlim[1] - y0b
				y1a <- ylim[2] + leftlim[1] - y1a
			}
			
			
			# Draw left branch
			# Rectangular branches
			if(type=="rect"){
				if(!horiz){
					segments(x0 = x0a, y0 = y0a, x1 = x0b, y1 = y0a)
					segments(x0 = x0b, y0 = y0a, x1 = x0b, y1 = y0b)
				}else{
					segments(y0 = x0a, x0 = y0a, y1 = x0b, x1 = y0a)
					segments(y0 = x0b, x0 = y0a, y1 = x0b, x1 = y0b)
				}
			# Triangular branches
			}else{
				if(!horiz){
					segments(x0 = x0a, y0 = y0a, x1 = x0b, y1 = y0b)
				}else{
					segments(y0 = x0a, x0 = y0a, y1 = x0b, x1 = y0b)				
				}
			}
			# Draw right branch
			# Rectangular branches
			if(type=="rect"){
				if(!horiz){
					segments(x0 = x0a, y0 = y0a, x1 = x1a, y1 = y0a)
					segments(x0 = x1a, y0 = y0a, x1 = x1a, y1 = y1a)
				}else{
					segments(y0 = x0a, x0 = y0a, y1 = x1a, x1 = y0a)
					segments(y0 = x1a, x0 = y0a, y1 = x1a, x1 = y1a)
				}
			# Triangular branches
			}else{
				if(!horiz){
					segments(x0 = x0a, y0 = y0a, x1 = x1a, y1 = y1a)
				}else{
					segments(y0 = x0a, x0 = y0a, y1 = x1a, x1 = y1a)
				}
			}
			
			# continue left branch
			plotdend(dend[[1]], type = type, horiz = horiz, new = F, leafpos = leafpos, xlim=xlim, ylim=ylim, xscale = xscale, yscale = yscale, left = left, revx=revx, revy=revy, normleafpos = F)
			# continue right branch
			plotdend(dend[[2]], type = type, horiz = horiz, new = F, leafpos = leafpos, xlim=xlim, ylim=ylim, xscale = xscale, yscale = yscale, left = left + attr(dend[[1]], "members"), revx=revx, revy=revy, normleafpos = F)
		}
	}

	# How reordering is to be done (rows)
	if(is.null(Rowv)){ # Default reordering
		Rowv <- rowMeans(x, na.rm = na.rm)
	}
	# Obtaining the dendrogram (rows)
	if(inherits(Rowv, "dendrogram")){ # Custom provided dendrogram
		dend.row <- Rowv
		order.row <- order.dendrogram(Rowv)
	}else if(!identical(Rowv, NA)){ # Default dendrogram
		# Dendrogram for each row
		dend.row <- as.dendrogram(hclustfun(distfun(x)))
		dend.row <- reorderfun(dend.row, Rowv)
		order.row <- order.dendrogram(dend.row)
	}else{ # Omit dendrogram (Rowv is NA)
		dend.row <- NULL
		order.row <- 1:nrow(x)
	}
	# Plot the dendrogram (cols)
	if(!is.null(dend.row) & !identical(Rowv, NA)){
		plotdend(dend.row, horiz=T, xlim=ylim, ylim=leftlim, new=F, type=type, revy=T)
	}

	# How reordering is to be done (cols)
	if(is.null(Colv)){ # Default reordering
		Colv <- colMeans(x, na.rm = na.rm)
	}
	# Obtaining the dendrogram (cols)
	if(inherits(Colv, "dendrogram")){ # Custom provided dendrogram
		dend.col <- Colv
		order.col <- order.dendrogram(Colv)
	}else if(!identical(Colv, NA)){ # Default dendrogram
		# Dendrogram for each col
		dend.col <- as.dendrogram(hclustfun(distfun(t(x))))
		dend.col <- reorderfun(dend.col, Colv)
		order.col <- order.dendrogram(dend.col)
	}else{ # Omit dendrogram (Colv is NA)
		dend.col <- NULL
		order.col <- 1:ncol(x)
	}
	# Plot the dendrogram (cols)
	if(!is.null(dend.col) & !identical(Colv, NA)){
		plotdend(dend.col, xlim=xlim, ylim=toplim, new=F, type=type)
	}
	
	# Scaling (copied as in stats::heatmap for compatibility, default value 'none' though)
	scale <- if (symm && missing(scale)) "none"
	else match.arg(scale)
	if (scale == "row") {
		x <- sweep(x, 1L, rowMeans(x, na.rm = na.rm), check.margin = FALSE)
		sx <- apply(x, 1L, sd, na.rm = na.rm)
		x <- sweep(x, 1L, sx, "/", check.margin = FALSE)
	}
	else if (scale == "column") {
		x <- sweep(x, 2L, colMeans(x, na.rm = na.rm), check.margin = FALSE)
		sx <- apply(x, 2L, sd, na.rm = na.rm)
		x <- sweep(x, 2L, sx, "/", check.margin = FALSE)
	}	
	
	# Finding color bins for the values according to the interval number in palette
	intervalmat = matrix(findInterval(x=x, vec=valseq), nrow=nr, ncol=nc)
	
	xmatseq = seq(from=xlim[1], to=xlim[2], length.out=nc+1)
	# Reversing the y-coordinates, otherwise the heatmap will be upside-down
	# y-axis runs [ymin, ymax] in different direction than rows are counted (first row is top-left corner)
	ymatseq = rev(seq(from=ylim[1], to=ylim[2], length.out=nr+1))
	
	# Reorder rows and columns for the color intervals according to the respective dendrograms
	intervalmat <- intervalmat[rev(order.row), order.col] # Notice that y-axis was reversed
	
	# Plotting the rectangles for the heatmap itself
	for(row in 1:nr){
		for(column in 1:nc){
			rect(xleft=xmatseq[column], ybottom=ymatseq[row], xright=xmatseq[column+1], ytop=ymatseq[row+1], col=col[intervalmat[row, column]], border=border[row, column], lty=lty[row, column], lwd=lwd[row, column])
		}
	}
	
	# Plot names for rows and columns
	# Custom row name vector provided
	if(!missing(namerows) & length(namerows)>1){
		rownames(x) <- namerows
		namerows <- T
	}
	# No column names, using a sequence of column indices
	else if(is.null(rownames(x))){
		rownames(x) <- 1:nrow(x)
	}
	if(is.null(rownames(x))) rownames(x) <- 1:nrow(x)
	rownam <- rownames(x)[order.row]
	xseq <- rep((rightlim[2] + rightlim[1])/2, times=(nrow(x)))
	yseq <- seq(from=ylim[1], to=ylim[2], length.out=(nrow(x)+1))
	yseq <- unlist(lapply(1:(length(yseq)-1), FUN=function(z) mean(c(yseq[z], yseq[z+1]))))
	rowtext <- list(xseq = xseq, yseq = yseq, rownam = rownam)
	# Plot names for rows
	if(namerows){
		textfun(xseq=xseq, yseq=yseq, labels=rownam, type="row")
	}
	# Custom column name vector provided
	if(!missing(namecols) & length(namecols)>1){
		colnames(x) <- namecols
		namecols <- T
	}
	# No column names, using a sequence of column indices
	else if(is.null(colnames(x))){
		colnames(x) <- 1:ncol(x)
	}
	colnam <- colnames(x)[order.col]
	xseq <- seq(from=xlim[1], to=xlim[2], length.out=(ncol(x)+1))
	xseq <- unlist(lapply(1:(length(xseq)-1), FUN=function(z) mean(c(xseq[z], xseq[z+1]))))
	yseq <- rep((bottomlim[2] + bottomlim[1])/2, times=(ncol(x)))
	coltext <- list(xseq = xseq, yseq = yseq, colnam = colnam)
	# Plot names for columns
	if(namecols){
		textfun(xseq=xseq, yseq=yseq, labels=colnam, type="col")
	}
	
	# Option values may be useful for plotting color key, annotations etc, so return invisibly
	invisible(list(
		x = x, 
		xmatseq = xmatseq,
		ymatseq = ymatseq,
		xlim = xlim,
		ylim = ylim,
		leftlim = leftlim,
		rightlim = rightlim,
		toplim = toplim,
		bottomlim = bottomlim,
		dend.row = dend.row, 
		dend.col = dend.col, 
		order.row = order.row, 
		order.col = order.col, 
		rowtext = rowtext, 
		coltext = coltext, 
		colors=col, 
		nbins = nbins, 
		valseq = valseq, 
		intervalmat = intervalmat
	))
}

# Function for plotting 'hmap' colour key
hmap.key <- function(
	# The invisible object returned by 'hmap' function for plotting parameters
	h,
	# Coordinates where the key box should be
	x0 = h$leftlim[1],
	x1 = h$leftlim[2],
	y0 = h$toplim[1],
	y1 = h$toplim[2],
	# Limits for the x-axis of the color key
	xlim = range(h$valseq),
	# y-axis ratio between value labels and the colors in the key box with top:bottom
	ratio = 0.5,
	# ratio in y-axis for value ticks
	tick = 0.1,
	# Values to annotate
	at = seq(from=min(h$valseq), to=max(h$valseq), length.out=5),
	# Should the key box be bounded ("o" is yes, "c" just colors box, "n" omits box)
	bty = "c",
	# Text size
	cex = 0.5,
	# Text position
	pos = 3#,

	# Additional parameters
	#...
){
	# If user wants to cut color key to specific range
	xseq <- seq(from=x0, to=x1, length.out=length(h$valseq[h$valseq >= xlim[1] & h$valseq <= xlim[2]]))
	h$colors <- h$colors[h$valseq >= xlim[1] & h$valseq <= xlim[2]]
	h$valseq <- h$valseq[h$valseq >= xlim[1] & h$valseq <= xlim[2]]
	# Compute y-axis middlepoints
	ymid <- y1-ratio*(y1-y0) # Middle point between color box and annotation box
	ytick <- y1-ratio*(y1-y0)-tick*(y1-y0) # Further tick annotated values from middle point downwards
	# Key colors
	for(i in 1:(length(xseq)-1)){
		rect(xleft=xseq[i], xright=xseq[i+1], ybottom=ymid, ytop=y1, col=h$colors[i], border=NA)
	}
	# Black boundary boxes
	if(bty=="o"){
		rect(xleft=x0, xright=x1, ybottom=y0, ytop=y1, border="black", col=NA)
	}else if(bty=="c"){
		rect(xleft=x0, xright=x1, ybottom=y1-ratio*(y1-y0), ytop=y1, border="black", col=NA)
	}
	# Value ticks and annotations
	where <- findInterval(at, vec=h$valseq)
	where <- where[!where==0]
	for(w in 1:length(where)){
		text(x=xseq[where[w]], y=y0, pos=pos, labels=round(at[w],3), cex=cex) # Text labels above lower text boundary box
		segments(x0=xseq[where[w]], x1=xseq[where[w]], y0=ymid, y1=ytick, col="black")
	}
	
	invisible(list(where = where, xseq = xseq, x0 = x0, x1 = x1, y0 = y0, y1 = y1, ymid = ymid, ytick = ytick))
}


# Function for annotating rows and columns in plots by 'hmap'
hmap.annotate <- function(
	# The invisible object returned by 'hmap' function for plotting parameters
	h,
	
	# ROWS
	# Groups per ROW to annotate, each unique value is assigned to a unique instance of the color palette (preferably a factor)
	rw,
	# Count of labels for ROWS
	rw.n = length(unique(rw)),
	# Color palette for ROWS
	rw.col = rainbow(rw.n, start=0.05, end=0.5),
	# If user desired plotted rectangles instead of pch-based symbols, use the widths and heights indicated here
	rw.wid,
	rw.hei,
	# Symbols for ROWS
	rw.pch,
	# x and y-axis coordinates to where the symbols ought to be plotted
	# ROWS
	rw.x = rep(min(h$rightlim), times=length(h$rowtext$xseq)),
	rw.y = h$rowtext$yseq,
	rw.shift = c(0.02, 0),
	# Position for a pre-generated row annotation legend, if NA or NULL then this legend is omitted
	#rw.pos = c(h$xlim[1],h$ylim[1]),

	# COLUMNS
	# Groups per COLUMN to annotate, each unique value is assigned to a unique instance of the color palette (preferably a factor)
	cl,
	# Count of labels for COLUMNSs
	cl.n = length(unique(cl)),
	# Color palette for COLUMNs
	cl.col = rainbow(cl.n, start=0.55, end=1.0),
	# If pch is present, the symbols in it are used for annotating the unique labels through recycling, if omitted then use pars rwid & rhei and clid & chei instead
	# If user desired plotted rectangles instead of pch-based symbols, use the widths and heights indicated here
	cl.wid,
	cl.hei,
	# Symbols for COLUMNS
	cl.pch,
	# x and y-axis coordinates to where the symbols ought to be plotted
	# COLUMNS
	cl.x = h$coltext$xseq,
	cl.y = rep(max(h$bottomlim), times=length(h$coltext$yseq)),
	cl.shift = c(0, -0.02),
	# Position for a pre-generated row annotation legend, if NA or NULL then this legend is omitted
	#cl.pos = c(h$xlim[1],h$ylim[1]),

	# Additional parameters
	...
){
	if(!missing(rw)){
		#if(missing(rw.wid) & missing(rw.hei)){
		if(!missing(rw.pch)){
			points(rw.x + rw.shift[1], rw.y + rw.shift[2], pch=rw.pch, col=rw.col[as.numeric(as.factor(rw))][h$order.row], ...)
		}else{
			for(i in 1:(length(h$ymatseq)-1)) {
				#rect(xleft=min(rw.wid), xright=max(rw.wid), ytop=h$ymatseq[i], ybottom=h$ymatseq[i+1], col=rw.col[as.numeric(as.factor(rw))][h$order.row][i], border=rw.col[as.numeric(as.factor(rw))][h$order.row][i], ...)
				rect(xleft=min(rw.wid), xright=max(rw.wid), ytop=h$ymatseq[i], ybottom=h$ymatseq[i+1], col=rw.col[as.numeric(as.factor(rw))][h$order.row][i], density=NA, ...)
			}
		}
		#if(!missing(rw.pos)){
		#
		#}
	}
	if(!missing(cl)){
		#if(missing(cl.wid) & missing(cl.hei)){
		if(!missing(cl.pch)){
			points(cl.x + cl.shift[1], cl.y + cl.shift[2], pch=cl.pch, col=cl.col[as.numeric(as.factor(cl))][h$order.col], ...)
		}else{
			for(i in 1:(length(h$xmatseq)-1)) {
				rect(xleft=h$xmatseq[i], xright=h$xmatseq[i+1], ytop=max(cl.hei), ybottom=min(cl.hei), col=cl.col[as.numeric(as.factor(cl))][h$order.col][i], density=NA, ...)
			}
		}
		#if(!missing(cl.pos)){
		#
		#}
	}

}


# Extend range -function (based on 'extendrange' from grDevices) with forced symmetry around a point
extendsymrange <- function(
	x, # extendrange x
	r = range(x, na.rm=T), # extendrange r
	f = 0.05, # extendrange f
	sym = 0 # Axis of symmetry
){
	ex <- extendrange(x = x, r = r, f = f)
	ran <- max(abs(ex - sym))
	c(sym - ran, sym + ran)
}

# Give in vector (or matrix/data.frame) of y-values, and obtain smart jittering for separating the measurements on the x-axis (obtain x-axis values)
smartjitter <- function(
	# Original data vector (or matrix or data.frame)
	x,
	# Quantile splits
	q = seq(from=0, to=1, length.out=10),
	# Type of jittering;
	# type == 1: consecutive jitters in same bin are values 0,0.1,0.2,0.3, ...
	# type != 1: consecutive jitters in the same bin are alternating -1 coefficiented values 0,-0.1,0.2,-0.3,0.4, ...
	type = 1,
	# Amount of jittering per overlapping values
	amount = 0.1,
	# Jittering function for values that reside in the same split bin
	jitterfuncs = list( 
		# type == 1 option function
		function(n){
			(1:n)/(1/amount)
		},
		# type != 1 option function
		function(n){
			(((-1)^c(0:(n-1)))*(0:(n-1)))/(1/amount)
		}
	),
	jits = jitterfuncs[[type]]
	)
{
	w <- as.numeric(cut(x, breaks=c(-Inf,unique(quantile(x, probs=q, na.rm=T)),Inf)))
	jit <- unlist(
		apply(
			do.call("rbind", 
				lapply(unique(w), FUN=function(z) 
					{ 
						res <- rep(NA, times=length(w)); 
						res[which(w==z)] <- jits(length(which(w==z))); 
						res 
					} 
				)), 
		MARGIN=2, FUN=function(y) 
			ifelse(all(is.na(y)), 0 , y[!is.na(y)]))
	)
	jit
}
