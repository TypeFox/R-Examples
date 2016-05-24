################################################################################
## PLOT.mvformula: Plot functions for mvabund objects                          #
## (multivariate abundance data)                                               #
################################################################################

default.plot.mvformula <- function(	x,
				y, 
				type="p", 
				main, 
				xlab,
  				ylab="Abundances", 
				col=if(type=="bx") "white" else "black", 
				fg="grey",
  				pch=1, 
				las=1, 
				write.plot="show", 
				filename="plot.mvabund",
  				n.vars=if(any(is.na(var.subset))) 12 else length(var.subset), 
				overall.main="",
  				data=parent.frame(), 
				var.subset=NA, 
				xvar.select=TRUE,
				n.xvars=NA,
#				n.xvars=if(any(is.na(xvar.subset))) min(3, sum(!is.interaction)) else length(xvar.subset),
				xvar.subset=NA, 
				scale.lab="ss", 	
				t.lab="t", 
				mfrow=NULL, 
        mfcol=NULL, 
				shift=TRUE, 
				border="black",
  				all.labels=FALSE, 
				keep.window=FALSE, 
				ask, ... ) {

#Kill Function if a BoxPlot has been passed
if (type=="bx") { stop("\nERROR: It is not ideal to use a boxplot for this data type\n") }

# If there is only one response variable this is changed to TRUE
# (overwritten independent of passed value)

dev <- dev.list()
dev.name <- getOption("device")

if(is.null(dev.name))
	stop("Make sure that the 'device' option has a valid value, e.g. 'options(device = 'windows')'. Allowed values here are 'windows', 'win.graph', 'x11', 'X11'.")

# if(!(any(dev.name == c("windows", "win.graph", "x11", "X11")) ) )
#   stop("Make sure that the 'device' option has a valid value, e.g. 'options(device = 'windows')'.
#   Allowed values here are 'windows', 'win.graph', 'x11', 'X11'.")

m <- match.call(expand.dots = FALSE)

dots  <- m$...

logWarn <- FALSE
if(!is.null(dots$log)){
	# dots$log <- NULL
	if(regexpr("y", dots$log ) !=-1) {
		dots$log <-  sub("y","", dots$log)
		logWarn <- TRUE
	}
}

usex <- TRUE
if (length(dots)>0){ 
	# Delete arguments in ... that are defined lateron and cannot be used twice in
	# the plot function.
	deactive <- c("axes", "cex.axis")  
	deactivate <- (1:length(dots))[names(dots) %in% deactive ]
	
	for (i in length(deactivate):1) {
		dots[ deactivate[i] ] <- NULL 			#fixed the [[]] compile problem, could reduace functionallity.
	}

	dots 	<- lapply(dots, eval, parent.frame())

	if(!is.null(dots$formula)) {
		x <- dots$formula 
		dots$formula <- NULL
		usex <- FALSE
	}

	subset.char <- dots$subset
} else subset.char <- NULL

miss.varsubset <- missing(var.subset)
# Change logical var.subset to numerical var.subset, if necessary. Note that
# NA values are logical as well, but should be excluded here.
if(!miss.varsubset){
	if(is.logical(var.subset) & any(!is.na(var.subset)))
	var.subset <- which(var.subset[!is.na(var.subset)])
}

miss.varsubset <- !is.numeric(var.subset) # If this function is called within
# another, the missing function could be tricked out.
var.subset  <- as.vector(var.subset)

miss.xvarsubset <- missing(xvar.subset)
if(!miss.xvarsubset){
	if(is.logical(xvar.subset) & any(!is.na(xvar.subset)))
	var.subset <- which(xvar.subset[!is.na(xvar.subset)])
}
miss.xvarsubset <- any(is.na(xvar.subset)) || any(is.null(xvar.subset))

if ( (missing(x) & usex) ) stop("'formula' is missing or incorrect")

if(!inherits(x, "formula") ) {
# if no formula is passed a formula of supposedly passed
# matrices/vectors/data.frames will be tried to build.
# If there is only one mvabund object, it is assumed to be the response.
# If there are two mvabund objects, x (the first argument passed) is the
# response, !note that this is different in plot.mvabund
# The response cannot be a data.frame, when both objects are data.frames an
# error will be produced.
# If none is an mvabund object, x is the response, if it is not a data.frame.
	if (missing(y)) { stop("A response and independent variables need to be provided for a formula.") }

	if( (is.mvabund(y) & ! is.mvabund(x)  ) | is.data.frame(x)) { 
		if(is.data.frame(y)) stop("a data.frame is not accepted as response for building the formula")
		targs <- match.call(call = sys.call(which = 1), expand.dots = FALSE)
		mvabund.formula <- x <- as.formula(paste(deparse(targs$y,width.cutoff = 500),"~",deparse(targs$x,width.cutoff = 500), sep=""))

	# That might result in troubles if either targs$y or targs$x has the name "x"
	# try to remove x (formula) earlier after giving it the name mvabund.formula
	} else { 
		targs <- match.call(call = sys.call(which = 1), expand.dots = FALSE)
		mvabund.formula <- x <- as.formula( paste(deparse(targs$x,width.cutoff = 500),"~", deparse(targs$y,width.cutoff = 500), sep="") )
	}
	# if(length(x) != 3)   stop("'formula' is incorrect")
} else mvabund.formula	<- x

resp <- attr(terms(mvabund.formula),"response")

term.labels <- attr(terms(mvabund.formula),"term.labels")

if(resp) {
	if(length(mvabund.formula) != 3) stop("'formula' is incorrect")
	mvabund.object.1<-as.matrix(model.response(model.frame( mvabund.formula, data=data)))
	nam <- paste(deparse( mvabund.formula[[2]]), ".", sep="")

	if(!is.null(colnames( mvabund.object.1 ))){
		if(any(is.na(suppressWarnings(as.numeric(sub( nam ,"", colnames( mvabund.object.1 ))))))){
			colnames( mvabund.object.1 ) <- sub( nam ,"", colnames( mvabund.object.1 ))
		}
	}
} else stop("'formula' has no response")

mvabund.object.1 <- mvabund(mvabund.object.1, neg=TRUE)
	
N <- NROW(mvabund.object.1)     # number of sites
p <- NCOL(mvabund.object.1)     # number of organism types

# Automatically correct n.vars if > p
n.vars <- min(n.vars,p)

if(length(term.labels)!=0){
	foo.ext <- extend.x.formula(mvabund.formula,return.interaction=TRUE, extend.term=TRUE)

	mvabund.formula  <- foo.ext$formula
	is.interaction   <- foo.ext$is.interaction
	# the length of this shows the number of x variables, including interaction
	# term and factors, but factors only appear once and not like in model.matrix
	# with (one column for every level) - 1
} else {
	is.interaction <- c()
}

# Exclude Interactions.
if(any(is.interaction)){

	mvabund.formula <- terms(mvabund.formula)[-which(is.interaction)]
	print(2)
	# terms.foo <- attr(model.frame(mvabund.formula), "terms")
	# facs <- attr(terms.foo,"factors")[,is.interaction, drop = FALSE]
	# datClass <- attr(terms.foo,"dataClasses")
	# for(faci in 1:NCOL(facs)){
	#     if(any(datClass[facs[,faci]>0] == "factor"))
	#       stop("Formulas with interactions that include factors cannot yet be plotted")
	# }
}

mf <- model.frame(mvabund.formula, data=data)
expl.data <- as.data.frame(mf[, -1])

# This is a matrix, ie factors are split to dummy variables
# not used for the actual plot, but for getting assign, etc
pExpl <- length(is.interaction) # NCOL(expl.data) + sum(is.interaction)
if(pExpl == 0) stop("The model is empty. There's nothing to plot.")

miss.xlab <- missing(xlab) 
if (miss.xlab) xlab <- NULL 

miss.main <- missing(main) 
if (miss.main)  {

	nam <- paste(deparse( mvabund.formula[[2]]), ".", sep="")
	# main <- colnames(as.matrix( model.frame( mvabund.formula)))[1:p]
	main <- colnames(mvabund.object.1)

	if(!is.null(main)){
		if(any(is.na(suppressWarnings(as.numeric( sub(nam,"", main) ))))){
			main <- colnames(mvabund.object.1) <- sub(nam,"", main)
		} else colnames(mvabund.object.1) <- main
	}

	# main <- colnames(mvabund.object.1) <- colnames(as.matrix(
	#   model.frame( mvabund.formula)))[1:p]
	n.mn <- nchar(main, type="chars")
	if(any(n.mn > 14))  main[n.mn > 14] <- paste(substr(main[n.mn > 14], 1, 7),"...\n...",
							substr(main[n.mn > 14], n.mn[n.mn > 14]-7, n.mn[n.mn > 14]) )
	
} else if (length(main)==p){
	colnames(mvabund.object.1) <- main

} else if (!all(is.na(var.subset)) & length(main)==length(var.subset)){
	colnames(mvabund.object.1)[var.subset]<- main
	main <- colnames(mvabund.object.1)

} else if (length(main)==1) {
	main <- rep(main,times =p)

} else stop("'main' must be a vector consisting of one or ", as.character(p), " characters" )

if (length(main)==1) { main <- rep(main,times= p) }

######### BEGIN edit var.subset, n.vars and mvabund.objects #########
# subset allows double variables
var.subset      <-  as.vector(var.subset)
var.subset.dim  <-  length(var.subset)

if (miss.varsubset) {
	if(!any(length(ylab)==c(1,p))) stop("the length of 'ylab' is not correct")
	sum.mvabund.object.1 <-
	t(mvabund.object.1) %*% matrix(1,ncol=1,nrow=N)

	# Find abundance ranks OF MVABUND.OBJECT.1.
	var.subset<-order(sum.mvabund.object.1, decreasing = TRUE)

	# Ensure no more than n.vars in var.subset.
	if (n.vars < length(var.subset)) var.subset <- var.subset[1:n.vars]

	var.subset.dim <- length(var.subset)

} else if ( p<max(var.subset) ) {
	# Some dimension checks for subset, n.vars.
	stop ("You have given an invalid 'var.subset'.")  

} else if (n.vars!=var.subset.dim) {
	n.vars <- var.subset.dim
}  

mvabund.object.1 <- mvabund.object.1[,var.subset, drop=FALSE] 

if(length(ylab)==p)  { 
	ylab <- ylab[var.subset] 
} else if (length(ylab)==1) ylab <- rep(ylab, times=n.vars)

main<- main[var.subset]
######### END edit var.subset, n.vars and mvabund.objects #########

mvabund.formulaUni <- formulaUnimva( mvabund.formula, var.subset=var.subset, split.x=TRUE, intercept=0 )

######### BEGIN edit xvar.subset, n.xvars etc #########

#DW changes, 30/10/14.  17/8: moved this outside of variable selection loop
#only define default n.xvars here, once is.interaction has been defined 
if(is.na(n.xvars))
{
  if(any(is.na(xvar.subset)))
    n.xvars = min(3, sum(!is.interaction))
  else
    n.xvars = length(xvar.subset)
}


if(xvar.select & pExpl>1) #begin varaible selection (if required)
{
	if(is.list(pch))
  {
		pch2 <- list()
		
		for(i in 1:length(xvar.subset))
    {
			pch2[[xvar.subset[i] ]] <- pch[[i]]
		}
		pch <- pch2
	}
	
	if(is.list(col))
  {
		col2 <- list()
	
		for(i in 1:length(xvar.subset))
    {
			col2[[xvar.subset[i]]] <- pch[[i]]
		}
		col<- col2
	}

  if(miss.xvarsubset)
  {
	  # DW change, 17/8/15: changed to a warning if no n.xvars passed. 
    if(any(names(dots)==n.xvars)==FALSE)
        warning(paste("No sQuote(n.xvars) passed so subset selection will be for",n.xvars,"variables"))
	
	   # Find the n.xvars independent variables with the highest average R^2.
	   xvar.subset <- best.r.sq(formula= mvabund.formula, # mvabund.formulaUni,
				var.subset=var.subset, n.xvars= min(n.xvars, pExpl))
     xvar.subset=xvar.subset$xs #to keep the indices only

		if (is.null(xvar.subset))
    {
		  xvar.subset <- 1:pExpl
    	warning("xvar.subset could not be found")
	  }
  }
  else if (length(xlab)==length(xvar.subset))
  {
    warning("The values of 'xlab' are taken to be corresponding to 'xvar.subset'")
    lab <- xlab
  	xlab <-  rep("",times=pExpl)
	 	xlab[xvar.subset] <- lab
		# It is required that xlab has a name for every column of x,
		# if xvar.subset is passed, xlab of the same length is accepted.
	}
  else xvar.subset <- 1:pExpl
}
else xvar.subset <- 1:pExpl

pExpl <- length(xvar.subset)
######### END edit xvar.subset, n.xvars etc #########

########### factor plot #################
is.all.factor <- TRUE
for(i in xvar.subset) {
	if (!is.factor(expl.data[,i])) is.all.factor <- FALSE
}

if(is.all.factor) {

	# If the indep. var are all factors, draw a factor plot.
	m[[1]]<- NULL
	m$...  <- NULL
	m$data <- NULL
	m$all.labels <- NULL
	m$x    <- NULL
	m$y    <- NULL
	m$xvar.select <- m$n.xvars <- m$xvar.subset <- m$var.subset <-
	m$nvars <- NULL

	m <- lapply(m, eval, parent.frame())
	lterm <- length(attr(terms(mvabund.formula),"term.labels"))

	if(is.null(dots$legendtitle) & ( lterm==1 | lterm == pExpl)) {
		m$legend.title <- attr(terms(mvabund.formula),"term.labels")[xvar.subset]
	} else if(!is.null(dots$title)){
		if(length(dots$title)>1) dots$title <- dots$title[xvar.subset]
	}

	if(shift) message("Overlapping points were shifted along the y-axis to make them visible.")


        cat("\n \tPIPING to 1st MVFACTOR PLOT \n") 
	   do.call( "default.plotMvaFactor", c(list(x=mvabund.object.1, y=expl.data[,xvar.subset, drop=FALSE],transformation="no"), m, dots))
  return(invisible())
}
########### END factor plot #################

# print a warning if the y-axis was wanted to be log-transformed - this should
# be done inside the formula
if(logWarn) warning("logarithmic axis is only possible on the x axis in 'plot.mvformula'")

########## whether to write the plot ##################################

if (write.plot!="show") {
	if (write.plot=="eps" | write.plot=="postscript")  {
		postscript(paste(filename,".eps", sep="") )
	} else if (write.plot=="pdf") {
		pdf(paste(filename,".pdf", sep="") )
	} else if (write.plot=="jpeg" ) {
		jpeg(paste(filename,".jpeg", sep="")) 
	} else if (write.plot=="bmp" ){
		bmp(paste(filename,".bmp", sep=""))
	}
	# Specify the window where to draw the plot.
	dev.curr <- dev.cur()
	on.exit( dev.off(which = dev.curr) )
}
###########

mvabund.formula <-  mvabund.formulaUni

###########

if (miss.xlab) xlab <- attr(terms(mvabund.formula[[1]]),"term.labels")
else
{ 
	llab <- length(attr(terms(mvabund.formula[[1]]),"term.labels"))

	# xlab needs to include names for every column of x.
	# not only the ones chosen with x.var.subset
	if (length(xlab)!= llab & length(xlab)>1)
		stop("'names.ind' must be a vector consisting of ",
			as.character(llab), " characters")
}  

############################################################################
# multiple plots are produced using the formula                            #
# according to R's defaults in the the univariate case                     #
############################################################################

if ((dim(expl.data)[1])!=N ) {stop ("You have not given independent data with appropriate dimension.")}

######### BEGIN establish row and column sizes #########


if (is.null(mfcol))
{
  #DW changes, 17-8-15: define mfrow here now that n.xvars is defined, and make it length one if one variable to plot
  if(length(n.xvars[xvar.subset])==1)
    mfrow=n.vars
  else
    mfrow=c(min(5,n.vars), min(3, n.xvars[xvar.select]))
  # end change
  perwindow <- mfrow
	mfr <- TRUE
} else {
	perwindow <- mfcol
	mfr <- FALSE
}

if(length(perwindow)==1) {

	if(any(c(pExpl, n.vars)==1)){
		columns <-ceiling(sqrt(perwindow))
		rows        <-columns-1
		if (rows*columns<perwindow) rows<-columns
	
	} else if (( (pExpl<7) & n.vars<5 & (pExpl*n.vars)<=perwindow ) | (pExpl>1 & perwindow%%pExpl == 0) ) {
		columns 	<- pExpl
		rows 		<- ceiling(perwindow/columns)
	
	} else {
		tim <- ceiling((pExpl*n.vars)/perwindow)
		if (pExpl > tim){
			rows 	<- ceiling(pExpl/tim)
			columns <- ceiling(perwindow/rows)
		
		} else {
			columns <- ceiling(n.vars/tim)
			rows 	<- ceiling(perwindow/columns)
		}
	}
} else {
	rows <- perwindow[1]
	columns	<- perwindow[2]
	perwindow <- rows*columns
}

if(any(c(pExpl, n.vars)==1)) {
	if (n.vars==1) all.labels <- TRUE
	
	if(missing(ask)) { 
		if (perwindow < pExpl*n.vars) { ask <- TRUE }
		else ask <- FALSE 
	}

} else {
	if(missing(ask)) { 
		if (n.vars>columns | pExpl>rows) { ask <- TRUE }
		else ask <- FALSE 
	}
}  

opar <- par( "mfrow", "mfcol", "mar", "ask") 

if((keep.window | all(opar$mfrow==c(rows,columns)))& write.plot=="show")
	opar$mfrow <- opar$mfcol<- NULL

if(write.plot=="show") on.exit(par(opar))

if(write.plot=="show" & is.null(dev)) {
	if (columns > (rows+1)){
		width   <- 16
		height  <- max(rows*width/columns,5)
	
	} else {
		height  <- 11
		width   <- max(height*columns/rows,4)
	}
	
	dev.off()
	do.call(dev.name, args=list(height=height,width=width))
}
######### END establish row and column sizes #########

if (all(par("mar")== c(5.1, 4.1, 4.1, 2.1)) & (rows*columns>1) )  {
	if (overall.main=="") side3 <- 1.7 else side3 <- 1
#	par(mar=c(3,4,side3,1)+0.1)
	par(mar=c(1.5,1.75,side3,1)+0.1)
}

#Set outer margin of the image (in text lines), and margin line on which labels are drawn.
#DW change, 4/5/15
#par(oma=c(1,1,2,1))
par(oma=c(2,2,2,1),mgp=c(1.75,0.5,0))


######### BEGIN plot #########
t.lab <- substr(t.lab,1,1)

if(!t.lab %in% c("o", "t")) stop("You have passed an invalid 't.lab'")

if (t.lab=="o") {
	# Obtain the transformation-axis ticklabels and tickmark values.
	# at the moment transAxis only works with the original name (so 'x' must be used)
	# and not if there are indices used, then the normal axis is drawn
	trAxis 	<- try(transAxis(x))
	if(class(trAxis)=="try-error"){
		axes <- TRUE
		yaxislabs <- yaxvalues <- NULL
	} else {

		if (trAxis$trans ) {  
			axes   <- FALSE 
			opp <- par("mgp")
			par(mgp=c(3,0.2,0))
			if(write.plot=="show") on.exit(par(opp), add =TRUE)

		} else axes <- TRUE

		yaxislabs <- trAxis$yaxlabs[var.subset]
		yaxvalues <- trAxis$yaxvalues[var.subset]
	}


} else if (t.lab=="t") {
	axes <- TRUE
	yaxislabs <- yaxvalues <- NULL

} else stop("undefined value in 't.lab'")


if (write.plot=="show") { aske <- ask  } 
else { aske <- FALSE }

if (!mfr) {
	rowsmfr <- rows
	colsmfr <- rows <- columns
	columns <- rowsmfr
}

if(any(c(pExpl, n.vars)==1)) { 
	layoutmat <- 1:min(pExpl*n.vars, perwindow)
} else {
	layoutmat <- (1:(min(n.vars,rows)*min(pExpl,columns)))
	layouthelp <- layouthelp2 <-matrix(0, ncol=columns,nrow=rows)
}

if(any(c(pExpl, n.vars)==1)) {

	layouthelp <- matrix( c(layoutmat,
	rep.int(0,times=max(c(0,rows*columns-max(layoutmat)))) ), ncol=columns,nrow=rows, byrow=TRUE )

	uselayout2 <- FALSE

} else {
	layouthelp[1:min(n.vars,rows),1:min(pExpl,columns)] <-
		matrix(layoutmat, nrow=min(n.vars,rows),ncol=min(pExpl,columns), byrow=TRUE )

	if(pExpl>columns & (pExpl%%columns)!=0) {
		layouthelp2[(1:min(n.vars,rows)),1:(pExpl%%columns)] <-
				matrix((1:(min(n.vars,rows)*(pExpl%%columns))),
					nrow=min(n.vars,rows),ncol=(pExpl%%columns), byrow=TRUE )
		uselayout2 <- TRUE
	} else uselayout2 <- FALSE
}


if (!mfr) {
	layouthelp <- t(layouthelp)
	if (uselayout2)  layouthelp2 <- t(layouthelp2)
}


if( write.plot!= "show")  {
	dev.set(which = dev.curr)
	# Specify the window where to draw the plot.
}


if (! all(is.null(c(mfcol, mfrow))) | is.null(dev) | write.plot!="show") {
	layout(layouthelp)
	nonewwind <- FALSE
	curr.layout <- layouthelp
} else { 
	nonewwind <- all.labels <- TRUE
}   


winxlevel <- ceiling(pExpl/columns)
winylevel <- ceiling(n.vars/rows)
	
if(pExpl==1) winylevel 	<-1
if(n.vars==1) winxlevel <-1
displab1x <- (!all.labels & pExpl==1)

if(any(scale.lab == "ss")){

	mxAll <- sapply( 1:n.vars, function(yi)	max(mvabund.object.1[,yi],na.rm=TRUE))
	wh.max <- which( abs( mxAll - max(mxAll, na.rm=TRUE)) < 1e-06  )[1]

	mnAll <- sapply( 1:n.vars, function(yi)min(mvabund.object.1[,yi],na.rm=TRUE))
	wh.min <- which( abs( mnAll - min(mnAll, na.rm=TRUE)) < 1e-06  )[1]

        ylimleft <- mnAll[wh.min]-mnAll[wh.min]/10
        ylimright <- mxAll[wh.max]+mxAll[wh.max]/10
	if(is.null(yaxvalues[[wh.max]])) {
		if(write.plot=="show" & is.null(dev)) {
			do.call(dev.name, args=list(height=height,width=width))
		} else { 
			do.call(dev.name, args=list())
		}	
	
#		dev.set(which = dev.cur())
cat("\n \tPIPING TO 1st PLOT FORMULA FEATURE \n")
#DW changes, 4/5/15:making this initial plot invisible. Still shows as a blank page though!? Need to remove...
#		do.call( "plotFormulafeature", c(list(mvabund.formula[[wh.max]],data=data, axes=axes, type=type[1],ylim=c(0,mxAll[wh.max]+mxAll[wh.max]/10), xvar.subset=1,ask=FALSE), dots ))
		do.call( "plotFormulafeature", c(list(mvabund.formula[[wh.max]],data=data, axes=axes, type="n",xaxt="n",yaxt="n",bty="n",ylim=c(ylimleft,ylimright), xvar.subset=1), xlab="",ylab="", dots ))
#		mtext(overall.main, outer = TRUE, cex = 0.9*par("cex.main"), col=par("col.main"), line = 1) 

		yaxTic <- axTicks(2)
		yaxLab <- as.character(yaxTic)
#		dev.off()

	} else {
		yaxTic <- yaxvalues[[wh.max]]
		yaxLab <- yaxislabs[[wh.max]]
	}

	yaxvalues <- yaxislabs <- vector("list", length=n.vars)
	
	for(yi in 1:n.vars) {
		yaxvalues[[yi]] <- yaxTic
		yaxislabs[[yi]] <- yaxLab
	}
}

par(ask=aske)
# if( shift & type != "bx" )
#  message("Overlapping points were shifted along the y-axis to make them visible.")

for (yi in 1:winylevel) {

	if(any(c(pExpl, n.vars)==1) | nonewwind) {
		rowsj <- 1:n.vars
	
	} else {         
		rowsj <- (1:rows) + rows*(yi-1)
		rowsj <- rowsj[rowsj<=n.vars]
	}

	
	for (xi in 1:winxlevel) {

		if (!nonewwind) {
			if (xi==winxlevel & uselayout2) {
				layout(layouthelp2) 
				curr.layout <- layouthelp2
			} else { 
				layout(layouthelp)
				curr.layout <- layouthelp
			}
		}
	
		if(any(c(pExpl, n.vars)==1) | nonewwind) {
			colsj <- xvar.subset
		} else {
			colsj <- (1:columns) + columns*(xi-1)
			colsj <- colsj[colsj<=pExpl]  
			colsj <- xvar.subset[colsj] 
		}

		# Get max value for y axis.
		max.y <- max(max(mvabund.object.1,na.rm=TRUE))
		min.y <- min(min(mvabund.object.1,na.rm=TRUE))

		for(yj in rowsj) {


			if(mfr) {
				if( all.labels | yj==rowsj[length(rowsj)] | (pExpl==1 & yj>n.vars-columns)) {
					xlabyj <- xlab # [colsj]
				} else  { 
					xlabyj <- ""
				}
 
				if(displab1x & (yj%%columns !=1)) ylabj <- "" 
				else ylabj <- ylab[yj]
	
			} else {
				xlabyj <- xlab # [colsj]
				if( all.labels | yj==rowsj[1] | (pExpl==1 & yj>n.vars-columns)) {
					ylabj <- ylab[yj] 
				} else ylabj <- ""
			}

			if(any(scale.lab == "ss")) {
		#		ylimVal <- c(0,mxAll[wh.max]+mxAll[wh.max]/10)
				ylimVal <- c(min.y, max.y)
			} else {
                #		mx <- max(mvabund.object.1[,yj],na.rm=TRUE)	
	        #		ylimVal <- c(0,mx+mx/10)
                                ylimVal <- c(min.y, max.y)
			}

cat("\n \tPIPING TO 2nd PLOT FORMULA FEATURE \n")
			do.call( "plotFormulafeature", c(list(mvabund.formula[[yj]], data=data,	ylab=ylabj, xlab=xlabyj, ask=aske, yaxis.ticks=yaxvalues[[yj]],	yaxis.labs=yaxislabs[[yj]], fg=fg, cex.xaxis=par("cex.axis"), cex.yaxis=par("cex.axis"), axes=axes , col=col, pch=pch, type=type, las=las, main= main[yj], ylim=ylimVal, border=border, xvar.subset=colsj, cex.lab=par("cex.lab"),all.labels= if(mfr) all.labels else TRUE), dots ))

      if (yj==rowsj[1] & overall.main!=""  & par("oma")[3]>=1) {
				mtext(overall.main, outer = TRUE, cex = 0.9*par("cex.main"), col=par("col.main"), line = 1) 
			}
		}
	}
}

if(xvar.select == TRUE &  pExpl < NCOL(expl.data) & miss.xvarsubset ) {
	# if(miss.xvarsubset)
	tmp <- " \n(the variables with the highest average R^2)" # else {
	#   tmp <- " (user selected)" }
	tmp2 <- paste(xlab[xvar.subset],  collapse = ", ") 
	message("Only the independent variables ", tmp2 , " were included in the plot(s)", tmp, "." )
}

if(n.vars < p & miss.varsubset) {
	#   if(miss.varsubset)
	tmp <- " \n(the variables with highest total abundance)"  # else {
	#    tmp <- " (user selected)"  }
	tmp2 <- paste(colnames(mvabund.object.1), collapse = ", ") 
	message("Only the response variables ", tmp2 , " were included in the plot(s)", tmp, ".")
}

if(!is.null(subset.char))  message("Only the subset ", subset.char,
	" of the cases was included in the plot(s) (user selected).")

pmfg <- par("mfg")
######### END plot #########

if(keep.window & !(pmfg[1]==rows & pmfg[2]==columns)) {

	if( (rows-pmfg[1]+columns-pmfg[2])==1) {
		# ie if only one free plot left.
		par(mfrow=c(rows,columns),mfg=c(rows,columns))
		return(invisible())
	} else if(columns==1) { 
		par(mfrow=c(rows,columns),mfg=c(pmfg[1]+1,1))
		return(invisible())
	} else if(rows==1) {
		par(mfrow=c(rows,columns),mfg=c(1, pmfg[2]+1))
		return(invisible())
	}

	if(max(curr.layout)== curr.layout[pmfg[1], pmfg[2]]) {
		# ie he current.layout is completely used.
		if(curr.layout[2,1]-curr.layout[1,1]==1) {
			wh.next <- which(curr.layout==0)
			mfgrows <- wh.next[1]%%rows
			if(mfgrows==0) mfgrows <- rows
			par(mfrow=c(rows,columns),mfg=c(mfgrows,ceiling(wh.next[1]/rows)))
		
		} else {
			wh.next <- which(curr.layout==0)
			mfgrows <- wh.next[1]%%rows
			if(mfgrows==0)  mfgrows <- rows
			par(mfcol=c(rows,columns),mfg=c(mfgrows,ceiling(wh.next[1]/rows)))
		}
	} else {
		# ie the curr layout has still at least one free place
		plot(0, axes=FALSE, type="n", xlab="", ylab="")
		# Check which is the next free place.
		pmfg2 <- par("mfg")
		if(pmfg2[1]>pmfg[1] | ( pmfg2[1]!= pmfg[1] & pmfg2[1]==1  ) ) {
			par(mfrow=c(rows,columns),mfg=c(pmfg2[1],pmfg2[2]))
		} else  par(mfcol=c(rows,columns),mfg=c(pmfg2[1],pmfg2[2]))
	}
	
}

}


