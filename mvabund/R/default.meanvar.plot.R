################################################################################
# PLOT.MEANVAR: Constructs mean/variance plots for an mvabund objects          #
# (multivariate abundance data)                                                #
################################################################################
default.meanvar.plot.mvabund <- function(x, 
				type="p", 
				log="xy", 
				main="", 
				sub=NULL,
  				xlab="",
				ylab="", 
				col, 
				pch=par("pch"), 
				asp=if(log=="xy") 1 else NA,
  				n.vars = NULL, 
				var.subset=NULL, 
				subset=NULL, 
				write.plot="show",
  				filename="plot.meanvar", 
				table=FALSE, 
				mfrow=1, 
				mfcol=NULL, 
				na.rm=TRUE,
  				add.trendline=FALSE, 
				trendline.col= 1, ... ) {

dev <- dev.list()
opp <- list() # get all the graphical parameters

dev.name <- getOption("device")

if(is.null(dev.name)) stop("Make sure that the 'device' option has a valid value,
	e.g. 'options(device = 'windows')'. Allowed values here are 'windows', 'win.graph', 'x11', 'X11'.")

#   if(!(any(dev.name == c("windows", "win.graph", "x11", "X11")) ))
#     stop("Make sure that the 'device' option has a valid value,
#     e.g. 'options(device = 'windows')'. Allowed values here are 'windows', 'win.graph', 'x11', 'X11'.")

if(all(is.null(c(mfrow, mfcol))))  op <- par("mfcol" ,"mfrow")

opp$new <- FALSE

if (any(is.na(x))) { 
	if(!na.rm ) { 
		stop("missing values in x")
	} else message("missing values in x were removed for calculations and plots")
}
	
mvabund.object  <- as.data.frame(x)
name.mvabund    <- (deparse(substitute(x)))

if(!is.null(subset)) {
	mvabund.object <- mvabund.object[c(subset),]
}

miss.varsubset <- missing(var.subset) | is.null(var.subset)
# Change logical var.subset to numerical var.subset, if necessary. Note that
# NA values are logical as well, but should be excluded here.
if(!miss.varsubset){
	if(is.logical(var.subset) & any(!is.na(var.subset)))
	var.subset <- which(var.subset[!is.na(var.subset)])
}
miss.varsubset <- !is.numeric(var.subset) 
# another, the missing function could be tricked out.
var.subset  <- as.vector(var.subset)

if (missing(n.vars) | is.null(n.vars)) {
	if(is.null(var.subset)) n.vars <- ncol(mvabund.object) 
	else n.vars <- length(var.subset)
}

rm(x)

# Assignments available for all the options.
N <- dim(mvabund.object)[1]     # number of sites
p <- dim(mvabund.object)[2]     # number of organism types
n.vars <- min(n.vars,p)   # n.vars is automatically corrected if >p


## BEGIN edit var.subset, n.vars and mvabund.objects
var.subset.dim <- length(var.subset)

if (miss.varsubset)  {
	# Arranging data to plot requested var.subset (default - n.vars most abund).
	if (n.vars>p) stop("You have passed an invalid number of variables to be included in the plot.")

	sum.mvabund.object <- t(mvabund.object)%*% matrix(1,nrow=N,ncol=1)
	
	# Find abundance ranks of the mvabund.object.
	var.subset <- order(sum.mvabund.object, decreasing = TRUE)

	# Ensure no more than n.vars var.subset.
	if (n.vars<length(var.subset)) var.subset<-var.subset[1:n.vars]

	var.subset.dim<-length(var.subset)

} else if (p<max(var.subset)) {
	stop ("You have passed an invalid var.subset") 

} else if (n.vars!=var.subset.dim) { 
	n.vars<- var.subset.dim
}

# Print the plot to the chosen medium.
if (write.plot!="show"){
	if (write.plot=="eps" | write.plot=="postscript") {
		postscript(paste(filename,".eps", sep="") )
	} else if (write.plot=="pdf") {
			pdf(paste(filename,".pdf", sep="") )  
	} else if (write.plot=="jpeg" ){
			jpeg(paste(filename,".jpeg", sep=""))    
	} else if (write.plot=="bmp" ){
			bmp(paste(filename,".bmp", sep=""))    
	} else if (write.plot=="png" ){
		png(paste(filename,".png", sep=""))
	}     
	on.exit( dev.off() )
} 

mvabund.object <- mvabund.object[,var.subset, drop=FALSE]
if (missing(col)) col <- "black"
if(length(col)==p)  col <- col[var.subset]
if(length(pch)==p)  pch <- pch[var.subset]
## END edit var.subset, n.vars and mvabund.objects

## BEGIN calculate overall means and variances for the abundances
mean.mvabund.overall <- colMeans(mvabund.object, na.rm=na.rm)
#mean.mvabund.overall <- mean(as.matrix(mvabund.object), na.rm=na.rm)
var.mvabund.overall  <- diag(var(mvabund.object, na.rm=na.rm))
## END calculate means and variances


## BEGIN calculate best fit 
if(add.trendline) {
	# Ensure that there are no NAs.
	nas <- is.na(mean.mvabund.overall+var.mvabund.overall)
	if(!all(nas)) {
		meanMVAnnas <- mean.mvabund.overall[!nas]
		varMVARnnas <- var.mvabund.overall[!nas]

		lm.fit <- lm(varMVARnnas ~ meanMVAnnas -1)

		#NO LONGER REQUIRED
		#beta.hat <- mean((meanMVAnnas-mean(meanMVAnnas))*
		#			(varMVARnnas-mean(varMVARnnas)))/mean((meanMVAnnas-mean(meanMVAnnas))^2)
		#NO LONGER REQUIRED
		#alpha.hat <- mean(varMVARnnas)-beta.hat*mean(meanMVAnnas)
	}
} else nas <- TRUE
## END calculate best fit 


## BEGIN edit label for x and y axis & cut away zeros for log-scale.
if ( regexpr("x",log)==-1 ) {

	meanMVA <- mean.mvabund.overall
	varMVAR <- var.mvabund.overall
	xlabel  <- "Mean [linear scale]"
} else {
	xlabel <- "Mean [log scale]"
	if (sum(mean.mvabund.overall<=0)>0) {
		message(sum(mean.mvabund.overall<=0),
		" mean values were <=0 and could not be included in the log-plot")
	}

	meanMVA <- mean.mvabund.overall[mean.mvabund.overall>0]
	varMVAR <- var.mvabund.overall[mean.mvabund.overall>0]

	if(length(col)==var.subset.dim)  col <- col[mean.mvabund.overall>0]
	if(length(pch)==var.subset.dim)  pch <- pch[mean.mvabund.overall>0]
}

if ( regexpr("y",log)== -1  ) {
	ylabel<-"Variance [linear scale]"  
} else {
	if (sum(var.mvabund.overall<=0)>0) {
		message(sum(var.mvabund.overall <= 0),
				" variance values were <=0 and could not be included in the log-plot")
	}

	ylabel  <- "Variance [log scale]"
	lvarMVAR <- length(varMVAR)
	if(length(col)==lvarMVAR)  col <- col[varMVAR>0]
	if(length(pch)==lvarMVAR)  pch <- pch[varMVAR>0]
	meanMVA <- meanMVA[varMVAR>0]
	varMVAR <- varMVAR[varMVAR>0]
}
## END edit label for x and y axis

if(is.null(main)) main <- paste(name.mvabund, "\n Mean-Variance Plot")
if(is.null(xlab)) xlab <- xlabel
if(is.null(ylab)) ylab <- ylabel

# Set up the window design.
if(!is.null(mfcol)) {
	mfrow <- mfcol
	mfr <- FALSE
}else if (!is.null(mfrow)) {
	mfr <- TRUE 
}else mfr <- NULL

if (length(mfrow)==1) {
	columns	<- ceiling(sqrt(mfrow))
	row <- columns-1
	
	if (row*columns<mfrow) row <- columns
	mfrow <- c(row,columns)
}

if(is.null(mfr)){
	row 	<- par("mfrow")[1]
	columns <- par("mfrow")[2]
}

if(write.plot=="show" & is.null(dev)) {
	
	if (columns > row){
		width 	<- 8
		height 	<- max(row*width/columns,5)
	} else {
		height 	<- 8
		width 	<- max(height*columns/row,4)
	}

	if(!is.na(asp)){
		if(regexpr("x", log)!= -1) r.mean <- log(range(meanMVA[meanMVA>0])) 
		else r.mean <- range(meanMVA)
	
		if(regexpr("y", log)!= -1) r.var <- log(range(varMVAR[varMVAR])) 
		else r.var <- range(varMVAR)
	
		hfac <- asp*(r.var[2]-r.var[1])/(r.mean[2]-r.mean[1])
	} else hfac <- 1

	if(!is.finite(hfac)) hfac <- 1
	if(hfac==0) hfac <- 1

	if(is.null(mfr)) dev.off()
			do.call(dev.name, args=list(height=height*hfac,width=width))
}
		
# Upon exiting the function, reset all graphical parameters to its value at
# the beginning, if the plot is drawn to a file, this will happen a anyway.
if( write.plot=="show" ) on.exit( par(opp) )

if(!is.null(mfr)){
	if (mfr) par(mfrow=mfrow)
	else par(mfcol=mfrow)
}

plot (meanMVA,varMVAR, xlab=xlab, ylab=ylab, type=type, main=main, sub=sub,log=log,col=col, pch=pch, asp=asp, ...)

###### START ADD TRENDLINE #######
if(add.trendline & !all(nas) ) {
	if(length(meanMVAnnas[!nas])>1) {

		if(regexpr("x", log)!= -1 | regexpr("y", log) != -1) {
cat("Lines 1\n")
			#mini <- max( min(meanMVAnnas[!nas]),1,(1-alpha.hat)/beta.hat )
			#maxi <- max( max(meanMVAnnas[!nas]),1,(1-alpha.hat)/beta.hat )
			#xlin <- exp(seq(log(mini),log(maxi), length.out = 30))
			#ylin <- alpha.hat +  xlin* beta.hat
			
			mini <- max( min(log(meanMVAnnas[meanMVAnnas > 0])) )
			maxi <- max( max(log(meanMVAnnas)) )
			xlin <- c(mini,maxi)
			lm.log <- lm(log(varMVAR[varMVAR > 0]) ~ log(meanMVA[meanMVA > 0]))
			ylin <- lm.log$coeff[1] +  xlin*lm.log$coeff[2]

			lines(exp(xlin),exp(ylin),col=trendline.col)
		} else {
cat("Lines 2\n")
			mini <- min(0, meanMVAnnas[!nas])
			maxi <- max(meanMVAnnas[!nas])
			beta.hat <- lm.fit$coeff
			lines(c(mini,maxi),c(mini*beta.hat,maxi*beta.hat),col=trendline.col)
		}
	}
}
###### END ADD TRENDLINE #######

if(n.vars < p) {
	if(miss.varsubset) tmp <- " \n(the variables with highest total abundance)" 
	else tmp <- " (user selected)"
	tmp2 <- paste(colnames(mvabund.object), collapse = ", ") 
	message("Only the response variables \n", tmp2, "\n", tmp ,
			" were included in the plot(s).")
}

if(!any(is.null(subset))) {
	message("Only the subset ", (deparse(substitute(subset))),
			" of the cases was included in the plot(s) (user selected).")
}

if (table) {
	MeanVarTable<-cbind(mean.mvabund.overall, var.mvabund.overall)
	colnames(MeanVarTable)<-c("Mean","Variance")
	return(MeanVarTable)
}

}
######## END MEANVAR.MVABUND FUNCTION #########

################################################################################
# PLOT.MEANVAR: Constructs mean/variance plots for an mvabund objects          #
# (multivariate abundance data)                                                #
# Separate means and variances are calculated for different samples as defined #
# by the formula, groups are chosen by factors in the independent data         #
# if there are no factors an overall Mean-Variance Plot is drawn               #
################################################################################
default.meanvar.plot.mvformula <- function(  x,  
					type="p", 
					log="xy", 
					main="",
  					sub=NULL, 
					xlab="", 
					ylab="", 
					col="", 
					pch=par("pch"),
  					asp=if(log=="xy") 1 else NA,
  					n.vars = NULL,
  					var.subset=NULL, 
					subset=NULL, 
					write.plot="show", 
					filename="plot.meanvar",
  					table=FALSE, 
					mfrow = if(index==0) c(1,1) 
						else if(overlay | write.plot=="show") 4 
						else c(min(5,n.vars),min(4,index)), 
					mfcol=NULL, 
					na.rm=TRUE,
  					add.trendline=FALSE, 
					trendline.col= 1, 
					cex.rel=FALSE,
  					overall.main=NULL, 
					overlay=TRUE, 
					all.labels=FALSE,
  					legend=FALSE, 
					legend.horiz, ... ) {

dev 	<- dev.list()
dev.name <- getOption("device")

if(is.null(dev.name)) stop("Make sure that the 'device' option has a valid value,
	e.g. 'options(device = 'windows')'. Allowed values here are 'windows', 'win.graph', 'x11', 'X11'.")

#   if(!(any(dev.name == c("windows", "win.graph", "x11", "X11")) ) )
#     stop("Make sure that the 'device' option has a valid value,
#     e.g. 'options(device = 'windows')'. Allowed values here are 'windows', 'win.graph', 'x11', 'X11'.")

dots <- match.call(expand.dots = FALSE)$...

formula.mvabund <- x
term.labels <- attr(terms(formula.mvabund), "term.labels")

if( attr(terms(formula.mvabund),"response")){
	mvabund.object <- as.data.frame(model.response(model.frame(formula.mvabund)) )
} else stop("formula has no response")

if(any(is.na(mvabund.object))) {
	if(!na.rm ){ stop("mvabund object contains NAs")
	} else message("missing values in the mvabund object were removed columnwise for calculations and plots")
}

miss.varsubset <- missing(var.subset) | is.null(var.subset)
# Change logical var.subset to numerical var.subset, if necessary. Note that
# NA values are logical as well, but should be excluded here.
if(!miss.varsubset){
    if(is.logical(var.subset) & any(!is.na(var.subset)))
    var.subset <- which(var.subset[!is.na(var.subset)])
}
miss.varsubset <- !is.numeric(var.subset) # If this function is called within
# another, the missing function could be tricked out.
var.subset  <- as.vector(var.subset)

rm(x)

N <- dim(mvabund.object)[1]     # number of sites
p <- dim(mvabund.object)[2]     # number of organism types

mvabund.colnames <- colnames(mvabund.object)
if (is.null(mvabund.colnames)) mvabund.colnames <- paste("Variable", 1:p)

variables <- eval( attr(terms(formula.mvabund),"term.labels") )
varlength <- length(variables)

nofactor <- TRUE
miss.n.vars <- missing(n.vars)

if (varlength==0) { # stop("formula has no explanatory variables")
	pExpl <- 0
} else {
	dataExpl <- as.data.frame(eval(parse(text=variables[1])) )
	pExpl <- NCOL(dataExpl)

	if( NROW(dataExpl)!=N ) stop("the dimensions of the variables do not match")

	if( varlength > 1 ) {
		for(i in 2:varlength) {
			dati <- eval(parse(text=variables[i]))
			if( NROW(dati)!=N ) stop("the dimensions of the variables do not match")
			dataExpl <- cbind(dataExpl, dati)
		}
	}

	ifa <- sapply(1:ncol(dataExpl), function(i) is.factor(dataExpl[,i]) )
	if(any(ifa)) nofactor <- FALSE
}

if(miss.n.vars){
	if(nofactor) { 
		nvars <- p
		if(any(is.null(var.subset))) n.vars <- p 
		else n.vars <- length(var.subset)
	} else {
		if(any(is.null(var.subset)) & overlay) { 
			nvars <- 12
		} else if (any(is.null(var.subset))) {  
			nvars <- 6
		} else nvars <- length(var.subset)
	}
}

n.vars <- min(n.vars,p)         # n.vars is automatically corrected if > p

if(!is.null(subset)) {
	# if there is only one factor and the col or pch has the length of the levels,
	# adjust the levels to new levels after subsetting.
	# ifa <- sapply(1:ncol(dataExpl), function(i) is.factor(dataExpl[,i]) )
	if (!nofactor) {
		if(sum(ifa)==1) {
			lev <- levels(factor(dataExpl[ ,which(ifa)]))
		}
		dataExpl <- dataExpl[subset,, drop = FALSE]

		if(sum(ifa)==1 & length(col)>=length(lev)) {
			lev2 <- levels(factor(dataExpl[ ,which(ifa)]))
			col <- col[which(lev %in% lev2)]
		}

		if(sum(ifa)==1 & length(pch)>=length(lev)) {
			lev2 <- levels(dataExpl[ ,which(ifa)])
			pch <- pch[which(lev %in% lev2)]
		}
	}

	mvabund.object <- mvabund.object[subset,, drop = FALSE]
	N <- dim(mvabund.object)[1]
}	


### BEGIN edit var.subset, n.vars and mvabund.objects
var.subset.dim <- length(var.subset)

if (miss.varsubset) {
	# Arrange data to plot requested var.subset (default - n.vars most abund).
	sum.mvabund.object <- t(mvabund.object) %*% matrix(1,nrow=N,ncol=1)
	
	if(na.rm & any(is.na(sum.mvabund.object))) {
		which <- which(is.na(sum.mvabund.object))
		for(i in which) sum.mvabund.object[i] <- sum(mvabund.object[,i], na.rm=TRUE)
	}

	# Find abundance ranks of mvabund.object.
	var.subset <- order(sum.mvabund.object, decreasing = TRUE)

	# Ensure no more than n.vars var.subset.
	if(n.vars < length(var.subset))  var.subset <- var.subset[1:n.vars]
	var.subset.dim <- length(var.subset)
	
} else if ( p < max(var.subset) ) {
	stop ("You have given an invalid var.subset")

} else if ( n.vars!=var.subset.dim )   {
	n.vars <- var.subset.dim
}


mvabund.object <- mvabund.object[ ,var.subset, drop=FALSE ]
mvabund.colnames <- mvabund.colnames[ var.subset ]
## END edit var.subset, n.vars and mvabund.objects


## BEGIN calculate means and variances for the abundances
index  <- max.length <- 0  # list-index

if(!nofactor) {

	expl.data.char <- attr(terms(formula.mvabund),"term.labels") 

	if (length(expl.data.char)!= pExpl) expl.data.char <- 1:length(pExpl)

	expl.cat.index <- levels <- mean.mvabund <- var.mvabund <- cex	<- iref <- n <- nlevels	<- list()

	for (i in 1:pExpl) {

		if(is.factor(dataExpl[,i]) ) {
			dataExpl[,i] <- factor(dataExpl[,i]) # necessary to get the updated Levels
	
			if (!na.rm& any(is.na(dataExpl[,i]))) {
				stop("independent variables object contain NAs")
			}

			index <- index+1  # the number of grouping variables
			nlevels[[index]] <- nlevels(dataExpl[,i])  
			levels[[index]]	<- levels(dataExpl[,i])

			max.length <- max(max.length, nlevels[[index]])
	
			expl.cat.index[[index]]	<-i

			mean.mvabund[[index]] <- matrix( ncol=n.vars, nrow=nlevels[[index]] )
			var.mvabund[[index]] <- matrix( ncol=n.vars, nrow=nlevels[[index]] )
			
			n[[index]] <- rep.int( NA,times=nlevels[[index]] )
			irind <- unclass( dataExpl[,i] )

			iref[[index]] <- irind[!is.na(irind)]
	
			for (j in 1:(nlevels[[index]]) ) { 
				var.mvabund[[index]][j,] <- diag(var(mvabund.object[iref[[index]]==j,,
										drop = FALSE], na.rm=na.rm))

#				mean.mvabund[[index]][j,] <- mean(as.matrix(mvabund.object)[iref[[index]]==j,,
				mean.mvabund[[index]][j,] <- colMeans(mvabund.object[iref[[index]]==j,, drop = FALSE], na.rm=na.rm)

				n[[index]][j] <- sum(iref[[index]]==j)
			}

			if (cex.rel) {
				cex[[index]] <- (log(n[[index]])+1)*2/max(log(n[[index]])+1)
			} else cex[[index]] <- 1
		}
	}
	expl.cat.index <- unlist(expl.cat.index)
}
## END calculate means and variances.

###### START 	#######	
if(!is.null(mfcol)) {
	mfrow <- mfcol
	mfr <- FALSE
} else if (!is.null(mfrow)) {
	mfr <- TRUE
} else mfr <- NULL
		
palet <- palette()

if (missing(col)) {
	
	if(max.length>15) {
		palette(rainbow(max.length))
		colr <- palette(rainbow(max.length))
	} else if (max.length>1) {

#		colr <- c("red", "darkgreen", "orange", "plum", "darkred", "darkblue",
#				"purple","rosybrown", "black", "green", "hotpink", "gold", "brown",
#					"lightblue","darkgrey")
		colr <- c(1:9)

		colr <- colr[1:max.length]

	} else colr <- "black"

	if (max.length > length(colr))	colr<- rep(colr, length.out=max.length)

} else if (max.length > length(col)) {
	message("As the suggested colors are not sufficient, they are used repeatedly.") 
	colr <- rep(col, length.out=max.length)
} else {
	colr <- col
}
col <- colr
	
if(length(pch)==1) { 
	pch <- rep(pch, times=max(max.length,1))
} else if (length(pch) < max.length) {
	message("As the symbols are not sufficient, they are used repeatedly.") 
	pch <- rep(pch, length.out=max(max.length,1))
}
	

## BEGIN edit label for x and y axis
if ( regexpr("x",log)==-1) {
	xmin	<-0
	xlabel	<-"[linear scale]"
} else {  
	xlabel	<-"[log scale]"
	xmin	<- NULL
}


if ( regexpr("y",log)==-1  ) {
	ymin	<-0
	ylabel	<-"[linear scale]"
} else {
	ylabel	<-"[log scale]"
	ymin	<- NULL
}
## END edit label for x and y axis

# BEGIN plot
if (write.plot!="show"){
	if (write.plot=="eps" | write.plot=="postscript") { 
		postscript(paste(filename,".eps", sep="") )  
	} else if (write.plot=="pdf") {
		pdf(paste(filename,".pdf", sep="") ) 
	} else if (write.plot=="jpeg" ) {
		jpeg(paste(filename,".jpeg", sep=""))  
	} else if (write.plot=="bmp" ) {
		bmp(paste(filename,".bmp", sep=""))    
	} else if (write.plot=="png" ){
		png(paste(filename,".png", sep=""))
	}        
	on.exit(dev.off())
	on.exit(palet, add=TRUE )
}

#mar <- par("mar")
namesub	<- deparse(terms(formula.mvabund)[[2]])


if (index==0) {
cat("START SECTION 1\n")
	## BEGIN calculate overall means and variances for the abundances
	# mean.mvabund.overall <- mean( as.matrix(mvabund.object), na.rm=na.rm )
	mean.mvabund.overall <- mean( mvabund.object, na.rm=na.rm )
	var.mvabund.overall  <- diag( var( mvabund.object, na.rm=na.rm ) )
	## END calculate overallmeans and variances   
	
	if (length(mfrow)==1)  {
		columns	<- ceiling(sqrt(mfrow))
		row <- columns-1
		if (row*columns<mfrow) row <- columns
		mfrow 	<- c(row,columns)
	}

	# Get all the graphical parameters
	opp <- par("ask", "mfrow", "mfcol")
#opp <- par("ask", "mar", "mfrow", "mfcol")

	if (is.null(mfr)) {
		# tmp       <- par("mfrow")
		row  <- par("mfrow")[1]
		columns	<- par("mfrow")[2]
		opp$mfrow <- opp$mfcol <- NULL
	}

	if(!is.null(mfr)){
		if (mfr) par(mfrow=mfrow)
		else par(mfcol=mfrow)
		row <- mfrow[1]
		columns	<- mfrow[2]
	}
	
	if (write.plot=="show") { 
		on.exit(par(opp))
		# Upon exiting the function, reset all graphical parameters to its
		# value at the beginning.
		on.exit(palet, add=TRUE )
	}

	if(write.plot=="show" & is.null(dev)) {
		if (columns > row){
			width 	<- 16
			height 	<- max(row*width/columns*1.2,5)
		} else {
			height 	<- 11
			width 	<- max(height*columns/row*0.83,4)
		}
		dev.off()
			
		if(!is.na(asp)){
			if(regexpr("x", log)!=-1){
				r.mean <- log(range(mean.mvabund.overall[mean.mvabund.overall>0]))
			} else r.mean <- range(mean.mvabund.overall)
	
			if(regexpr("y", log)!=-1){
				r.var <- log(range(var.mvabund.overall[var.mvabund.overall>0]))
			} else r.var <- range(var.mvabund.overall)

			hfac <- asp*(r.var[2]-r.var[1])/(r.mean[2]-r.mean[1])
		} else hfac <- 1

		if(!is.finite(hfac)) hfac <- 1
		if(hfac==0) hfac <- 1
	
		do.call(dev.name, args=list(height=height*hfac,width=width))
	}

	if ( regexpr("x",log)==-1  )  {
		meanMVA <- mean.mvabund.overall
		varMVAR <- var.mvabund.overall
	} else {
		if (sum(mean.mvabund.overall<=0,na.rm=na.rm)>0) {
			message(sum(mean.mvabund.overall<=0,na.rm=na.rm),
				" mean values were <=0 and could not be included in the log-plot")
		}
		meanMVA <- mean.mvabund.overall[mean.mvabund.overall>0]
		varMVAR <- var.mvabund.overall[mean.mvabund.overall>0]
	}

	if ((regexpr("y",log)!=-1 )& (sum(var.mvabund.overall<=0,na.rm=na.rm)>0)) {
		message(sum(var.mvabund.overall<=0,na.rm=na.rm),
				" variance values were <=0 and could not be included in the log-plot")
		meanMVA <- meanMVA[varMVAR>0]
		varMVAR <- varMVAR[varMVAR>0]
	}

#FIX TRENDLINE#
	if(add.trendline){
		# Ensure that there are no NAs.
		nas <- is.na(meanMVA+varMVAR)
		if(!all(nas)){
			meanMVAnnas <- meanMVA[!nas]
			varMVARnnas <- varMVAR[!nas]
			beta.hat <- mean((meanMVAnnas-mean(meanMVAnnas))*
						(varMVARnnas-mean(varMVARnnas)))/mean((meanMVAnnas-mean(meanMVAnnas))^2)
			alpha.hat <- mean(varMVARnnas)-beta.hat*mean(meanMVAnnas)
		}
	} else nas <- TRUE

	if (missing(main)) main <- "Overall mean-variance relation of abundances"
	if (missing(xlab)) xlab <- paste("Mean",namesub, xlabel)
	if (missing(ylab)) ylab <- paste("Variance",namesub, ylabel)

	plot( meanMVA, varMVAR,log=log, xlab=xlab, ylab=ylab,
			type=type,main=main,sub=sub, pch=pch, asp=asp, col=col, ... )

#FIXTRENDLINE##
	if(add.trendline & !all(nas)) {
		if(length(meanMVAnnas)>1) {
			if(regexpr("x", log)!= -1 | regexpr("y", log) != -1) {
				mini <- max( min(meanMVAnnas[!nas]),1,(1-alpha.hat)/beta.hat )
				maxi <- max( max(meanMVAnnas[!nas]),1,(1-alpha.hat)/beta.hat )
				xlin <- exp(seq(log(mini),log(maxi), length.out = 30))
				ylin <- alpha.hat +  xlin* beta.hat
				lines(xlin,ylin,col=trendline.col)
			} else {
				mini <- min(meanMVAnnas[!nas])
				maxi <- max(meanMVAnnas[!nas])
				lines(c(mini,maxi),c(alpha.hat+mini*beta.hat,alpha.hat+maxi*beta.hat),col=trendline.col)
			} 
		}
	}
	
	if (table) {
		MeanVarTable<-cbind(mean.mvabund.overall, var.mvabund.overall)
		colnames( MeanVarTable )<-c( "Mean", "Variance" )
	}
cat("FINISHED SECTION 1 \n")
} else {
cat("START SECTION 2 \n")
	opp <- par("ask", "mfrow", "mfcol") 
#opp <- par("ask", "mar", "mfrow", "mfcol") 

	if (length(mfrow)==1) { perwindow <- mfrow } 
	else {
		if(is.null(mfr)) { 
			mfr <-TRUE 				
			mfrow <- par("mfrow")
			if (! overlay)	opp$mfrow <- opp$mfcol <- NULL
		}
		rows <- mfrow[1]
		cols <- mfrow[2]
		perwindow <- rows*cols
	}
	
	if (write.plot=="show") {
		on.exit(par(opp))
		on.exit(palet, add=TRUE )
	}

	if (missing(xlab)) xlab <- paste("Mean", xlabel)
	if (missing(ylab)) ylab <- paste("Variance", ylabel)

	if (table) {

		MeanVarTable<-list()
		MVT.cn <- c("n", c(paste("Mean",  mvabund.colnames), paste("Var",
					mvabund.colnames))[rep(1:n.vars, each=2)+rep(c(0,n.vars),times=n.vars)])

		for (ind in 1:index) {
			if (length(n[[ind]])>1)
				MeanVarTable[[ind]]<-cbind(n[[ind]],(cbind(mean.mvabund[[ind]],
						var.mvabund[[ind]]))[,rep(1:n.vars, each=2)+rep(c(0,n.vars),
							times=n.vars)])
			else {
				MeanVarTable[[ind]]<-c(n[[ind]],(c(mean.mvabund[[ind]],
						var.mvabund[[ind]]))[rep(1:n.vars, each=2)+rep(c(0,n.vars),
							times=n.vars)])
				MeanVarTable[[ind]]<-t(as.matrix(MeanVarTable[[ind]]))
			}
			
			colnames(MeanVarTable[[ind]])<-MVT.cn
			rownames(MeanVarTable[[ind]])<- levels[[ind]] 
			names(MeanVarTable)[ind]<-expl.data.char[expl.cat.index[ind] ]
		}	
	}

	if (overlay) {
cat("Plotting if overlay is TRUE\n")
#		mar[c(3,4)] <- c(3.1,2)

		if (!missing(main)) {
			if(length(main)==1) main <- rep(main, times=index)
			else if (length(main)>0 & length(main)!=pExpl) stop("main must have length ", pExpl)
			else main <- main[expl.cat.index]
		}

		windows <- ceiling(index/perwindow)
		if (windows > 1 & write.plot=="show")  par(ask=TRUE)
	
		perwind<-min(index,perwindow)

		if (length(mfrow)==1) {

			if(perwind<=5) {
				cols <- perwind
				rows <- 1
			} else {
				rows <- ceiling(sqrt(perwind))
				cols <- rows-1
				if (cols*rows<perwind) rows <- cols
			}
		}
		
		if (legend){
			cols <- cols		
			# 2*cols to have place for the legend next to the plot.
			perwindleg <- perwind
		} else perwindleg <- perwind
		
		if(write.plot=="show" & is.null(dev)) {
			if (cols > rows){
				width 	<- 8
				height 	<- max(rows*width/cols,5)
			} else {
				height 	<- 8
				width 	<- max( height*cols/rows,4)
			}
			try(dev.off(), silent=TRUE)
			do.call(dev.name, args=list(height=height,width=width))
		}
		
		layoutmat <- c(1:(perwindleg), rep.int(0, times=(cols*rows)-(perwindleg)))
		layoutmat <- matrix(layoutmat , ncol=cols, nrow=rows, byrow=mfr)

#Checked - Set the window size for plot (need to check when multiple plots per window)

		
		for (l in 1:windows) {
			kchoice <- ((l-1)*perwind+1):((l-1)*perwind+perwind)
			kchoice <- kchoice[kchoice<(perwind+1)]
			ncoll <- rep.int(0,times=length(kchoice))
			
			for (ind in kchoice) ncoll[ind] <- min(ceiling(nlevels[[ind]]/(5+5*rows)),4)
			
			widths <- rep (c(4, max(ncoll)+1), length.out=cols)

#Check to see if commenting out this kills plot? Can't find relevance or point of using this?
#This bit of code is used to generate room for the legend. It creates a seperate plot for it.
			layout(layoutmat, widths= widths )
	
			for (ind in kchoice) { 
#				par(mar=mar)
				var.mvabund[[ind]][is.na(var.mvabund[[ind]])] <-0

				# Otherwise all the variances for only one group member are NA.
				mainlabel <- expl.data.char[expl.cat.index[ind] ]
				
				if (missing(main)) {
					main.ind <- paste("mean-var plot,", mainlabel, sep=" ")
				} else {
					main.ind <- mainlabel<- main[ind]
				}

				xmaxind<- max(mean.mvabund[[ind]], na.rm=na.rm)
				ymaxind<- max(var.mvabund[[ind]], na.rm=na.rm)

				if ( regexpr("x",log)==-1  ) {
					xminind <- min(mean.mvabund[[ind]], na.rm=na.rm) 
				} else  {
					xminind <- min(mean.mvabund[[ind]][mean.mvabund[[ind]]>0],na.rm=na.rm)
	
					if (sum(mean.mvabund[[ind]]<=0, na.rm=na.rm)>0) {
						message("using grouping variable ",mainlabel," ",
								sum(mean.mvabund[[ind]]<=0,na.rm=na.rm),
									" mean values were 0 and could 
										not be included in the log-plot")
					}
				}

				if ( regexpr("y",log)==-1  )  {
					yminind<- min(var.mvabund[[ind]], na.rm=na.rm)
				} else {
					yminind<- min(var.mvabund[[ind]][var.mvabund[[ind]]>0])
					
					if (sum(var.mvabund[[ind]]<=0,na.rm=na.rm)>0) {
						message("using grouping variable ",
								mainlabel," ", sum(var.mvabund[[ind]]<=0,na.rm=na.rm),
									" variance values were 0 and could not 
										be included in the log-plot")
					}
				}

				plot(c(xminind,xmaxind),c(yminind,ymaxind), type="n", asp=asp,
						log=log, xlab=xlab, ylab=ylab, main=main.ind,sub=sub, ... )
						
				j <- 1
				while (j< (nlevels[[ind]]+1)) {
					points(mean.mvabund[[ind]][j,],var.mvabund[[ind]][j,],
						col=colr[j],type=type,cex=cex[[ind]][1], pch=pch[j],...)

					j=j+1
				}
				
				if(add.trendline) {
					nas <- is.na(mean.mvabund[[ind]]+var.mvabund[[ind]])
					# Ensure that there are no NAs.
					meanindnnas 	<- c(mean.mvabund[[ind]])[!nas]
					varindnnas 	<- c(var.mvabund[[ind]])[!nas]
					
					if(is.numeric(meanindnnas+varindnnas) & length(meanindnnas)>1 ){
						# Check if nas consists only of FALSE values and more than one meanindnas.
						#beta.hat <- mean((meanindnnas-mean(meanindnnas))*
						#			(varindnnas-mean(varindnnas)))/
						#				mean((meanindnnas-mean(meanindnnas))^2)
						#alpha.hat <- mean(varindnnas)-beta.hat*mean(meanindnnas)
					
						if(regexpr("x", log)!= -1 | regexpr("y", log) != -1) {
							
							t.var <- log(varindnnas[varindnnas>0])
							t.mean <-log(meanindnnas[meanindnnas>0])

							mini <- min(t.mean)
							maxi <- max(t.mean)

							lm.fit <- lm(t.var~t.mean)
							beta.hat <- lm.fit$coeff
							lines(exp(c(mini,maxi)),exp(c(beta.hat[1]+mini*beta.hat[2],
											beta.hat[1]+maxi*beta.hat[2])),
												col=trendline.col)
						} else {
							lm.fit <- lm(varindnnas~meanindnnas -1)
							beta.hat <- lm.fit$coeff
							mini <- min(0,meanindnnas)
							maxi <- max(meanindnnas)
							lines(c(mini,maxi), c(mini*beta.hat, maxi*beta.hat),col=trendline.col)
						}				
					}
				}
				
				if (legend) {
#browser()
					#par(mar=c(0,0,mar[3],0))
				
					#plot(0,0,type="n",axes=FALSE,xlab="", ylab="")
					if(missing(legend.horiz)) {
						ncoll <- ceiling(nlevels[[ind]]/(5+5*rows))
					} else  ncoll <- 1
					
					if(ncoll==1) legloc <- "bottomright" 
					else legloc <- "topleft"
					
					cexl<-0.95
					
					if (ncoll>4) {
						ncoll<-4
						cexl <- min(0.9, 26/nlevels[[ind]])
					}
					
					leg <- levels[[ind]]
					
					tmp   <- suppressWarnings(as.numeric( leg ))
					natmp <- is.na(tmp)
				
					leg[natmp]  <- substr( leg[natmp], 1, (8/ncoll)+1 )
					leg[!natmp] <- zapsmall( as.numeric(leg[!natmp]), digits=(8/ncoll)+1 )

					# if(all(is.numeric(suppressWarnings(as.numeric( leg ))))) {
					# leg <- zapsmall(as.numeric(leg), digits=(8/ncoll)+1) } else {
					# leg <- substr(leg, 1,(8/ncoll)+1) }
					
					if(!is.null(dots$title)){
						legend( x=legloc, legend=leg, col=as.vector(colr[1:nlevels[[ind]] ]),
							ncol=ncoll, cex=cexl, pch=as.vector(pch[1:nlevels[[ind]] ]), ...)
					} else {
						if(length(term.labels)==1) title <- term.labels 
						else title <- NULL
						legend( x=legloc, legend=leg, col=as.vector(colr[1:nlevels[[ind]] ]),
							ncol=ncoll, cex=cexl, pch=as.vector(pch[1:nlevels[[ind]] ]),
							title = title, ...)
					}
				}
			}

			if (par("oma")[3]>=1)  {
				mtext(overall.main, outer = TRUE, cex = 1.1*par("cex.main"),col=par("col.main"))
			}
		}
	} else { 
cat("Plotting if overlay is FALSE\n")
#		if(!all.labels) mar[4] <- 0.5

		if (!missing(main)) {
			if(is.list(main)){
				if(length(main[[1]])!=p & !is.null(main[[1]]))
					stop("'main[[1]]' must be a character vector of length ",p)
				sublabel <- (main[[1]])[var.subset]
		
				if(!is.null(main[[2]]) & length(main[[2]])!=pExpl)
					stop("'main[[2]]' must be a character vector of length ",pExpl)
		
				mainlabel <- (main[[2]])[expl.cat.index]
			
			} else if(length(main)==1) {
				sublabel <- NULL
				mainlabel <- rep(main,times=n.vars)
			
			} else stop("'main' must be a character or a list of two character vectors containing ",
					"the names of the response variables and the names of the independent variables")
		} else {
			sublabel <- mvabund.colnames
			mainlabel <- expl.data.char[expl.cat.index ]	
		}	
	
		if (length(mfrow)==1) {
			if (index*n.vars <= mfrow) {
				cols <- index		#  index = the number of grouping variables
				rows <- ceiling(mfrow/cols) #  }
			} else {
				tim <- ceiling((index*n.vars)/mfrow)	# number of windows necessary
				rows <- ceiling(index/tim)
	
				if (index > tim) {
					cols <- ceiling(mfrow/rows)
				} else {
					if(legend) {
						rows <- 2 # One row will be required for the legend
						cols <- mfrow/2
					} else cols <- mfrow 	# rows = 1 when index <= tim
				}
			}
		}


		rowsleg <- rows
		if(legend) {
			rows <- rows-1 # else rowsleg <- rows
			legd <- 1
		} else legd <- 0
	
		if(!mfr) {
			colshelp <- cols
			rowshelp <- rowsleg
			cols <- rowshelp
			rowsleg <- colshelp
		}
	
		layoutmat <- matrix(0, ncol=cols, nrow=rowsleg) 
		timesrow <- rowsleg %/% (n.vars+ legd)
		timescol <- cols %/% index

		if(timesrow==0){
			layoutmat2 <- layoutmat
			submat <- rowsleg*(cols-(max(cols,index) %% index))
			layoutmat[,1:(cols-(max(cols,index) %% index)) ] <- 1:submat

			if(((n.vars+ legd)  %% rowsleg) ==0){
				submat2 <- ((n.vars+ legd)  %% rowsleg) *(cols-(max(cols,index) %% index))
				layoutmat2[((n.vars+ legd)  %% rowsleg),
						1:(cols-(max(cols,index) %% index)) ] <- 1:submat2
				uselayout2 <- TRUE
			} else uselayout2 <- FALSE
		} else {
			submat <- (n.vars+ legd)*(cols-(max(cols,index) %% index))
			uselayout2 <- FALSE

			for(i in 1:timesrow){
				layoutmat[(i-1)*(n.vars+ legd) + (1:(n.vars+ legd)),
						1:(cols-(max(cols,index) %% index)) ] <- (i-1)*(submat) + 1:submat
			}	
		}
	
		if(!mfr) {
			# mfcol was passed, change back now that layoutmat was calculated
			layoutmat <- t(layoutmat)
		
			if(timesrow==0){layoutmat2 <- t(layoutmat2) }
			cols <- colshelp
			rowsleg <- rowshelp
			# Avoid confusion and strange plots.
			all.labels <- TRUE
		}
		## END establish row and column sizes

		meanmax	<- meanmin <- varmax <- varmin <- rep.int(NA, times =index)

		for (ind in 1:index){
			var.mvabund[[ind]][is.na(var.mvabund[[ind]])] <-0
			# Otherwise all the variances for only one group member are NA.

			meanmax[ind]<- max(mean.mvabund[[ind]], na.rm=na.rm)
			varmax[ind]<- max(var.mvabund[[ind]], na.rm=na.rm)
			
			if ( regexpr("x",log)==-1  ) {
				meanmin[ind]<- min(mean.mvabund[[ind]], na.rm=na.rm)
			} else {
				meanmin[ind]<- min(mean.mvabund[[ind]][mean.mvabund[[ind]]>0])
				if (sum(mean.mvabund[[ind]]<=0,na.rm=na.rm)>0) {
					message("using grouping variable ",
							mainlabel," " , sum(mean.mvabund[[ind]]<=0,na.rm=na.rm),
								" mean values were 0 and could not be included in the log-plot")
				}
			}

			if ( regexpr("y",log)==-1 ) {
				varmin[ind]<- min(var.mvabund[[ind]], na.rm=na.rm)
			} else {
				varmin[ind]<- min(var.mvabund[[ind]][var.mvabund[[ind]]>0])
				if (sum(var.mvabund[[ind]]<=0,na.rm=na.rm)>0) {
					message("using grouping variable ",
							mainlabel," ", sum(var.mvabund[[ind]]<=0,na.rm=na.rm),
								" variance values were 0 and could not be included in the log-plot")
				}
			}
		}

		if(write.plot=="show" & is.null(dev)) {
			if (cols > rowsleg) {
				width 	<- 16
				height 	<- max( rowsleg*width/cols*1.2, 5)
			} else {
				height 	<- 11
				width 	<- max( height*cols/rowsleg*0.8,4)
			}
			try(dev.off(),silent=TRUE)
			
			if(!is.na(asp)) {
				m.min <- min(meanmin)
				m.max <- max(meanmax)
				v.min <- min(varmin)
				v.max <- max(varmax)
	
				if(regexpr("y", log)!= -1) {
					v.min <- log(v.min)
					v.max <- log(v.max)
				}
	
				if(regexpr("x", log)!=-1) {
					m.min <- log(m.min)
					m.max <- log(m.max)
				}

				hfac <- asp*(v.max-v.min)/(m.max-m.min)
				if(!is.finite(hfac)) hfac <- 1
				if(hfac==0) hfac <- 1
			} else hfac <- 1
			
			do.call(dev.name, args=list(height=height*hfac,width=width))
		}

		winxlevel <- ceiling(index/cols)
		winylevel <- ceiling(n.vars/rows)
		if ((winxlevel>1 | winylevel>1) & write.plot=="show") par(ask=TRUE)
		
		layout( layoutmat )

		for (yi in 1:winylevel) {
			rowsj <- (1:rows) + rows*(yi-1)
			rowsj <- rowsj[rowsj <= n.vars]
			
			if(uselayout2 & yi == winylevel) layout(layoutmat2)

			for (xi in 1:winxlevel) {
				
				colsj <- (1:cols) + cols*(xi-1)
				colsj <- colsj[colsj<=index]

				for (xj in colsj) {
						
					if (all.labels | xj == colsj[1]) ylabj <- ylab
					else ylabj <- "" 
						
					for (yj in rowsj){
#						if(!all.labels) {
#							if (yj==rowsj[1]) par(mar=c(3.1,mar[2:4]))
#							else par(mar=c(mar[1:2],2.6,mar[4]))

#						} else par(mar=mar)
							
						if (all.labels | yj == rowsj[length(rowsj)]) xlabj <- xlab
						else xlabj <- ""
						
						if (all.labels | yj == rowsj[1]) {
							main <- c( mainlabel[xj],"\n",sublabel[yj])
						} else main=c("\n",sublabel[yj])

						plot(c(meanmin[xj],meanmax[xj]),c( varmin[xj],varmax[xj] ),
								log=log, xlab=xlabj, ylab=ylabj, main=main, sub=sub,
									type="n", asp=asp,...)

						# Use the min and max of th matrices
						# to have the same plot dimensions for every x
						points(mean.mvabund[[xj]][,yj],var.mvabund[[xj]][,yj],
							col=as.vector(colr[1:nlevels[[xj]] ]) ,
							type=type,cex=cex[[xj]] , pch=as.vector(pch[1:nlevels[[xj]] ]),...)

						if(add.trendline) {
							nasj <- is.na(mean.mvabund[[xj]][,yj]+var.mvabund[[xj]][,yj])
							# Ensure that there are no NAs.
							meanindnnasj 	<- (mean.mvabund[[xj]][,yj])[!nasj]
							varindnnasj 	<- (var.mvabund[[xj]][,yj])[!nasj]

							if(is.numeric(meanindnnasj+varindnnasj) & length(meanindnnasj)>1) {
								# Check if nas consists only of FALSE values.
								beta.hat <- mean((meanindnnasj-mean(meanindnnasj))*
											(varindnnasj-mean(varindnnasj)))/
											mean((meanindnnasj-mean(meanindnnasj))^2)
								
								alpha.hat <- mean(varindnnasj)-beta.hat*mean(meanindnnasj)
								nas <- is.na(mean.mvabund[[xj]]+var.mvabund[[xj]])
								
								# Use the min and max of the matrices.
								meanindnnas <- (mean.mvabund[[xj]])[!nas]
								varindnnas <- (var.mvabund[[xj]])[!nas]
	
								if(regexpr("x", log)!= -1 | regexpr("y", log) != -1) {
									mini <- max( min(meanindnnas),1,(1-alpha.hat)/beta.hat )
									maxi <- max(max(meanindnnas),1,(1-alpha.hat)/beta.hat )
									xlin <- exp(seq(log(mini),log(maxi), length.out = 30))
									ylin <- alpha.hat +  xlin* beta.hat
									lines(xlin,ylin,col=trendline.col)
								} else {
									mini <- min(0, meanindnnas)
									maxi <- max(meanindnnas)
									lines(c(mini,maxi),c(alpha.hat+mini*beta.hat,
											alpha.hat+maxi*beta.hat),col=trendline.col)
								}
							}
						}
					}
						
					if (legend) {
						if(!mfr) {
#							par(mar=c(0,0.5,4.1,0.5)) 
							ncolmax <- 4
							ncolfact <- 9
						} else {
#							par(mar=c(0,4.1,0.5,0.5))
							ncolmax <- 3
							ncolfact <- 12
						}
	
						plot(0,0,type="n",axes=FALSE,xlab="", ylab="")
						cexl <- 0.9
						leg <- levels[[ind]]

						tmp   <- suppressWarnings(as.numeric( leg ))
						natmp <- is.na(tmp)
						leg[natmp]  <- substr( leg[natmp], 1, 5 )
						leg[!natmp] <- zapsmall( as.numeric(leg[!natmp]), digits=5 )

						if(missing(legend.horiz)) {
							ncoll <- ceiling(nlevels[[ind]]/min(ncolfact,24/rows))
							if (ncoll>ncolmax) {
								ncoll<-ncolmax
								cexl <- min(0.7, 26/nlevels[[ind]])
							}
							legend(title=mainlabel[xj], x="topleft",legend=leg,horiz=legend.horiz,
									col=as.vector(colr[1:nlevels[[xj]] ]), cex=cexl,ncol=ncoll,
										pch=as.vector(pch[1:nlevels[[xj]] ]), ... )	
						} else {
							ncoll <- 1

							if(length(leg)>27) {
								cexl <- 0.7

								tmp   <- suppressWarnings(as.numeric( leg ))
								natmp <- is.na(tmp)
								leg[natmp]  <- substr( leg[natmp], 1, 4 )
								leg[!natmp] <- zapsmall( as.numeric(leg[!natmp]), digits=4 )

								maxlev <- c(13, 26, 39)
							} else maxlev <- c(9,18,27)
	
							nlevind  <- length(leg)
							nlevind1 <- min(nlevind,maxlev[1])
							legend( title=mainlabel[xj], x="topleft",legend=leg[1:nlevind1],
							horiz=legend.horiz, col=as.vector(colr[1:nlevels[[xj]] ])[1:nlevind1],
							cex=cexl, ncol=ncoll,pch=as.vector(pch[1:nlevels[[xj]] ])[1:nlevind1], ... )
	
							if(nlevind > maxlev[1]) {
								nlevind1 <- min(nlevind,maxlev[2])
								legend( title=mainlabel[xj], x="left",
								legend=leg[(maxlev[1]+1):nlevind1], horiz=legend.horiz,
								col=as.vector(colr[1:nlevels[[xj]] ])[(maxlev[1]+1):nlevind1],
								cex=cexl, ncol=ncoll,
								pch=as.vector(pch[1:nlevels[[xj]] ])[(maxlev[1]+1):nlevind1], ... )
							}
	
							if(nlevind > maxlev[2]) {
								nlevind1 <- min(nlevind,maxlev[3])
								legend( title=mainlabel[xj], x="bottomleft",
								legend=leg[(maxlev[2]+1):nlevind1], horiz=legend.horiz,
								col=as.vector(colr[1:nlevels[[xj]] ])[(maxlev[2]+1):nlevind1],
								cex=cexl, ncol=ncoll,
								pch=as.vector(pch[1:nlevels[[xj]] ])[(maxlev[2]+1):nlevind1], ... )
							}
						}
					}
				}
				
				if (par("oma")[3]>=1) {	
					mtext(overall.main, outer = TRUE, cex = 1.1*par("cex.main"),col= par("col.main") )
				}
			}
		}
	}
cat("FINISHED SECTION 2 \n")
}

if(pExpl > index & index!=0) {
	tmp <- paste(expl.data.char[expl.cat.index],collapse = ", ")
	message("Only the independent variables \n", tmp,
			"\n (the factors) were included in the plot(s)." )
}

if(n.vars < p) {
	if(miss.varsubset) tmp <- " \n(the variables with highest total abundance)"   
	else tmp <- " (user selected)"

	tmp2 <- paste(colnames(mvabund.object), collapse = ", ") 

	message("Only the response variables \n", tmp2 ,
			"\nwere included in the plot(s)", tmp, ".")
}

if(!any(is.null(subset))) {
	message("Only the subset ", (deparse(substitute(subset))),
			"  of the cases was included in the plot(s) (user selected).")
}

if (table) return(MeanVarTable)

}

