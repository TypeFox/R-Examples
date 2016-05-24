################################################################################
# PLOTMVAFACTOR                                                                #
################################################################################
default.plotMvaFactor <- function(	x, 
				y, 
				type="p", 
				main="Abundances", 
				xlab,
  				ylab="Species", 
				col, 
				fg= "grey", 
				pch,
  				las=1, 
				write.plot="show", 
				filename="plot.mvabund", 
				n.vars= min(12,NCOL(x)),
  				overall.main, 
				var.subset=NA, 
				subset=NA, 
				transformation="log", 
				scale.lab="s",
  				t.lab="o", 
				mfrow = min(4,pExpl), 
				mfcol=NULL, 
				shift=TRUE, 
				border,
  				keep.window=FALSE, 
				ask, 
				legend = TRUE, 
				legend.horiz=FALSE,
  				legend.title, 
                                xaxis.labels, ... ) {


dev <- dev.list()
dev.name <- getOption("device")
#if (!is.null(dev)) dev.off() # close previous window
    
if(is.null(dev.name))
     stop("Make sure that the 'device' option has a valid value, e.g. 'options(device = 'windows')'.	Allowed values here are 'windows', 'win.graph', 'x11', 'X11'.")

#     if(!(any(dev.name == c("windows", "win.graph", "x11", "X11")) ) )
#       stop("Make sure that the 'device' option has a valid value,
#       e.g. 'options(device = 'windows')'. Allowed values here are 'windows', 'win.graph', 'x11', 'X11'.")

allargs <- match.call(expand.dots = FALSE)

dots <- allargs$...

if(!is.null(dots$log)){
     dots$log <- NULL
     warning("argument 'log' not implemented in 'plotMvaFactor'")
}

if (missing(x)) { stop("The mvabund object 'x' is missing.") }
if (missing(y)) { stop("The factors 'y' are missing.") }

if (write.plot!="show") {
    if (write.plot=="eps" | write.plot=="postscript") {
       	postscript( paste(filename,".eps", sep=""))
    } else if (write.plot=="pdf") {
	pdf(paste(filename,".pdf", sep=""))
    } else if (write.plot=="jpeg" ) {
	jpeg(paste(filename,".jpeg", sep=""))
    } else if (write.plot=="bmp" ) {
	bmp(paste(filename,".bmp", sep=""))
    } else if (write.plot=="png" ) {
	png(paste(filename,".png", sep=""))
    }
    
    # Close current window before proceeding to ensure plots are drawn to the right device.
    dev.curr <- dev.cur()
    on.exit( dev.off(which = dev.curr) )
}

if (length(dots)>0) {
# Delete arguments in ... that are defined lateron and cannot be used
# twice in the plot function
   deactive <- c("xlim", "ylim", "axes", "horizontal", "names")
   deactivate <- (1:length(dots))[names(dots) %in% deactive ]
        
   for (i in length(deactivate):1) {
       	dots[ deactivate[i] ]<-NULL #fixed the [[]] compile problem, could reduce functionality
    }
}	
    dots <- lapply(dots, eval, parent.frame(), parent.frame())
        
    if ("cex.lab" %in% names(dots)) clab <- dots$cex.lab
    else clab <- 1.4
        
    if ("col.lab" %in% names(dots)) colab <- dots$col.lab
    else colab <- par("col.lab")
        
    if ("cex" %in% names(dots)) cex <- dots$cex
    else cex <- 1.5

    if ("lwd" %in% names(dots)) lwd <- dots$lwd
    else lwd  <- 2 

    if ("cex.axis" %in% names(dots)) caxis <- dots$cex.axis
    else caxis <- 1.5

    if ("legend.text" %in% names(dots)) leg <- dots$legend.text
    else leg <- NULL 


if(! (type %in% c("p", "bx", "n"))) stop("The 'type' ", type, " is not implemented for factors as independent variables.")

mvabund.object <- as.matrix(unabund(x))
expl.data <- as.data.frame(y)

# Initiate n.vars before deleting x.
n.vars <- n.vars
if(!any(is.na(subset))) {
   mvabund.object <- mvabund.object[c(subset),, drop=FALSE]
   expl.data <- as.data.frame(expl.data[subset,])
}

N <- nrow(mvabund.object)     # number of sites
p <- ncol(mvabund.object)     # number of organism types

# Automatically correct n.vars if > p
n.vars <- min(n.vars,p)

if (any(c(p,N)==0)) stop("The mvabund object 'x' has invalid dimensions.")

miss.varsubset <- missing(var.subset)
# Change logical var.subset to numerical var.subset, if necessary. Note that
# NA values are logical as well, but should be excluded here.
if(!miss.varsubset) {
    if(is.logical(var.subset) & any(!is.na(var.subset))) var.subset <- which(var.subset[!is.na(var.subset)])
}
miss.varsubset <- !is.numeric(var.subset) # If this function is called within
# another, the missing function could be tricked out.
var.subset <- as.vector(var.subset)

pExpl <-  ncol(expl.data)

if(nrow(expl.data) !=N)
      stop("the dimension of the factor/s 'y' doesn't match the dimension of the mvabund object")

tmp <- NULL

# Check if all columns in x are factors.
for(i in 1:pExpl) if (!is.factor(expl.data[,i])) tmp <- c(tmp, i)

# Non-factors are deleted.
if(!is.null(tmp)) expl.data <- expl.data[,-tmp]

pExpl.old <- pExpl
pExpl <-  NCOL(expl.data)

miss.pch     <- missing(pch)
miss.col     <- missing(col)
miss.border  <- missing(border)
miss.legtitle  <-  missing(legend.title)
 
#rm(x)
mvabund.colnames <- colnames(mvabund.object)

#if(miss.legtitle) {	
#    if(is.factor(y)) factornames <- deparse( substitute(y))[1]
#    else factornames <-  colnames(expl.data)    
#    if(is.null( factornames)) factornames <- paste("factor", 1:pExpl, sep="")
#} else {
#    factornames <- legend.title
#    if(length(legend.title)==1) {
#	if(pExpl>1) factornames <- paste(legend.title, 1:pExpl, sep="")
#      	else factornames <-legend.title
#    } else if(length(legend.title)==pExpl) { 
#	factornames <- legend.title
#    } else if(length(legend.title)==pExpl.old) { 
#	factornames <- legend.title[-tmp]
#    } else stop("'legend.title' must be a character vector of the length of the number of factors or of length 1")
#}
    
if(is.null(mvabund.colnames)) mvabund.colnames <- paste("Variable", 1:p)

opp <- par("ask", "col.lab", "cex.lab", "mar", "fg", "mgp" ,"mfcol" ,"mfrow")
    
if(keep.window & write.plot=="show") opp$mfrow <- opp$mfcol <- NULL

# Avoid that a new window is drawn with the next plot after exiting the function.
if(!is.null(mfcol)) mfrow <- mfcol else opp$mfrow <- opp$mfcol <- NULL

if (length(mfrow)==1)  {
       	columns <- ceiling(sqrt(mfrow))
       	row     <- columns-1

      	if (row*columns<mfrow) row<-columns

       	mfrow = c(row,columns)
         	
	if(write.plot=="show") {
	   if(all(opp$mfrow==c(row,columns))) opp$mfrow <- opp$mfcol <- NULL
        }
} else {
       	row <- mfrow[1]
       	columns <- mfrow[2]
}

# Upon exiting the function, all graphical parameters are reset to its
# values at the beginning.
if(write.plot=="show") on.exit(par(opp))

if((write.plot=="show" & is.null(dev) ) & (!is.null(mfrow))) {        
   if (mfrow[2] > mfrow[1]){
       width   <- 16
       height  <- max(mfrow[1]*width/mfrow[2],5)*1.4 # *0.7	
   } else if (mfrow[2] == mfrow[1]){       
       if (legend == TRUE) {
           width <- 6 
           height  <- max(mfrow[1]*width/mfrow[2],5)*1.5 #0,7 
       }	   
       else {
           width <- 6
           height  <- max(mfrow[1]*width/mfrow[2],6)*1.4 # *0.7 
       }	   
   } else {
       height  <- 11
       width   <- max(height*mfrow[2]/mfrow[1],4) *0.4 # *1.4
   }
   # Close the old window first.
   dev.off()
   do.call(dev.name, args=list(height=height,width=width))
   devflag <- TRUE
} else devflag <- FALSE

#Change made for multiple plots in the window.
if (!is.null(mfcol)) par(mfcol=mfcol)
else if (!is.null(mfrow)) par(mfrow=mfrow)

if(missing(ask))  {
    if(pExpl > prod(mfrow)  &  write.plot=="show") { ask = TRUE } 
    else ask <- FALSE 
}

######### BEGIN edit var.subset, n.vars and mvabund.objects  #########
var.subset.dim <- length(var.subset)
if (!is.numeric(var.subset)) {
     if (n.vars>p) stop("You have passed an invalid number of variables 'n.vars' to be included in the plot.")
     sum.mvabund.object <- t(mvabund.object)%*%matrix(1,ncol=1,nrow=N)       
      # Find abundance ranks of mvabund.object.
        # Do some dimension checks for the subset.
        var.subset <- order(sum.mvabund.object, decreasing = TRUE)
        if (n.vars < length(var.subset)) var.subset <- var.subset[1:n.vars]
        
	# Ensure no more than n.vars in var.subset.
        var.subset.dim <- length(var.subset)

        # Arrange data to plot requested var.subset
        # (default is to the n.vars most abundand variables).
} else if (p<max(var.subset)) {
         stop ("You have passed an invalid 'var.subset'")
} else if (n.vars!=var.subset.dim) { 
	n.vars<- var.subset.dim
}
  	
mvabund.object   <- mvabund.object[,var.subset, drop=FALSE]
mvabund.colnames <- mvabund.colnames[var.subset]

#Get Min & Max or Transformation of Axis
min0 <- x[,var.subset]
tick.min <- min(min0[min0!=0],na.rm=TRUE)
tick.max <- max(x[,var.subset],na.rm=TRUE)

######### END edit var.subset, n.vars and mvabund.objects #########

########## BEGIN check transformation, t.lab, scale.lab ###########
if(substr(transformation, 1,1)=="n" | transformation==""){
      	transformation <- "no"
} else if(substr(transformation, 1,1)=="l"){ 
	transformation <- "log"
} else if(substr(transformation, 1,1)=="s" & transformation=="sqrt"){
      	transformation <- "sqrt"
} else if(substr(transformation, 1,1)=="s" & transformation=="sqrt4"){
      	transformation <- "sqrt4"
} else stop("You have passed an invalid 'transformation'")
   
scale.lab <- substr(scale.lab,1,1)
if(!scale.lab %in% c("r", "s")) stop("You have passed an invalid 'scale.lab'")
t.lab <- substr(t.lab,1,1)
if(!t.lab %in% c("o", "t")) stop("You have passed an invalid 't.lab'")
########### END check transformation, t.lab, scale.lab ############

# Get the variable numbers, which the abundances are plotted against.
# Make sure the plot starts at the top of y-axis going downwards.
y.axis <- rep(n.vars:1 , each=N)

# Get the max value for y axis before any transformations.
mx <-  max(mvabund.object, na.rm=TRUE)

# Do a 's'tandard plot: 0 is included in the x axis.
if ( scale.lab=="s") {
	xlim <- c(-0.1,mx+mx/20)
# Do an 'r'-plot: R's default is used
} else if ( scale.lab=="r") xlim <- NULL  
else {
	stop ("undefined value for 'scale.lab'")
}

#Get Min & Max or Transformation of Axis
#tick.min <- min(mvabund.object[mvabund.object!=0],na.rm=TRUE)
#tick.max <- max(mvabund.object,na.rm=TRUE)

######### START transformation #########
if ((transformation!="no") ) {
    if (t.lab=="o" ) {
	if (scale.lab=="s") ylim <- c(1,mx+mx/10) else ylim <- NULL
	if(devflag) do.call(dev.name, args=list(height=height,width=width))  
	else do.call(dev.name, args=list())

	if(any(c(mvabund.object)!=0) ) {
            suppressWarnings( plot( c(mvabund.object),type="n",axes=FALSE,xlab="", ylab="",ylim=NULL,log="y"))
            axis3 <- axTicks(2)
	} else {
	    plot( mvabund.object,y.axis,type="n",axes=FALSE,xlab="", ylab="", xlim=xlim,ylim=c(0,n.vars+0.5))
            # Get the transformation-labels for the third axis.
	    axis3 <- axTicks(2)
	    lengthax3 <- length(axis3)
	    if (lengthax3 > 5) axis3 <- axis3[-(5:(lengthax3-1))]
	}
	dev.off()
	ax3lab <-  as.character(axis3)
	
	if (scale.lab=="s") {
	    ax3lab <- c("0", ax3lab )
	    axis3 <- c(0,axis3) 
	}
    }

    ######### BEGIN transformation #########
    # Transform data, if required.
    if (transformation=="log") {
	if (max(mvabund.object,na.rm=TRUE)==0)
            stop("The mvabund object 'x' only consists of zero-abundances")
        minNon0 <- min(mvabund.object[mvabund.object!=0],na.rm=TRUE)
	mvabund.object <- log(mvabund.object+minNon0)-log(minNon0)

	# plot title: 'Abundances [log(y/min+1) scale]'
	transf.lab <- expression(paste("Abundances  ", bgroup("[", paste(log, bgroup( "(",frac(y,min)+1,")"), " scale"),"]")))
            
	if (t.lab=="o" ) {axis3 <- log(axis3+minNon0)-log(minNon0)}
        tick.min <- log(tick.min)
	tick.max <- log(tick.max)
     } else if (transformation=="sqrt4") {
         mvabund.object <- (mvabund.object)^0.25
	 # Define plot title:  'Abundances [y^{0.25} scale]'.
	 transf.lab <-  expression(paste("Abundances  ", bgroup("[", paste(sqrt(y,4)," scale") ,"]")))

	 if (t.lab=="o") {axis3<-(axis3)^0.25 }
         tick.min <- tick.min^0.25
	 tick.max <- tick.max^0.25
    } else if (transformation=="sqrt") {
         mvabund.object <- sqrt(mvabund.object)
	 # Define plot title: 'Abundances [y^{0.5} scale]'.
	 transf.lab <-  expression(paste("Abundances  ", bgroup("[", paste(sqrt(y)," scale") ,"]")))
	 if (t.lab=="o") { axis3 <- sqrt(axis3)  }
         tick.min <- sqrt(tick.min)
	 tick.max <- sqrt(tick.max)
    }
} else transf.lab <- 'Abundances'
######### END transformation #########

if( write.plot!= "show")  {
	# Avoid confusion where to draw the plot, if not specified the plot could
      	# be drawn to a window which is already open at the beginning.
      	dev.set(which = dev.curr)
}

######### BEGIN some calculations for better axis scaling #########
#if (missing(xlab)) 
   xlab <- transf.lab

# Get minimum value for correction of posx.
minmva <- min(mvabund.object,na.rm=TRUE)

# Get maximum value for y axis after transformation.
mx     <- max(mvabund.object,na.rm=TRUE)

#if (scale.lab=="s" ) {
#	xlim <- c(-0.1,mx+mx/20)
#	if(mx==0){
#		mx <- 1
#		xlim = c(-0.1, 0.1)
#	} else if (mx<0.8) {
#		potl <- -floor(log(mx)/log(10) )
#		mx   <- ceiling(mx *(10^potl)) /(10^potl)
#	} else  if (mx>6) {
#		
#		if (mx<10) potl <- - ceiling(log(mx)/log(10) ) 
#		else { potl <- - floor(log(mx)/log(10) ) }
#		
#		mx <- ceiling(mx *(10^potl)) /(10^potl)
#	} else mx <- ceiling(mx)

#	seque <- seq(from=0, to=mx, by=mx/10)
#	
#	if (mx > 1000) {
#		potence <- 10^(floor(log(mx)/log(10)))
#		sequenc <- as.character( round(seque/potence,digits=2 ))
#		xlab <- paste( xlab, " in ", potence)
#		sequenc <- as.character(round(seque,digits = 2 ) )
#	} else sequenc <- as.character(round(seque,digits = 2 ) )
#} else {
#	seque <- NULL
#	sequenc <- NULL
#	labels <- NULL
#	xlim <- NULL
#}
######### END some calculations for better axis scaling #########

if (scale.lab=="s" ) {
	xlim <- c(-0.1,mx+mx/20)
	
#	tick <- axisTicks(transform=transformation, max=tick.max, min=tick.min, tran.lab=t.lab)
	tick <- c(transform=transformation, max=tick.max, min=tick.min, tran.lab=t.lab)
#	seque <- tick$x.tic 
#	sequenc <- tick$x.ticlab
        seque <- NULL
	sequenc <- NULL
} else {
	seque <- NULL
	sequenc <- NULL
	labels <- NULL
	xlim <- NULL
}

##############################
######### BEGIN plot #########
##############################
  
############  either draw a boxplot  #################
if (type=="bx")  {
	y.axisat  <-  y.axis
	# variable numbers, the abundances are plotted against them
	# yaxis data should start at the top of y-axis going downwards, use at for this
	y.axis <- rep(1:n.vars , each=N)
	mvabnam <- substring(mvabund.colnames,first=1, last= max(12, nchar(mvabund.colnames, type="char") )  )

	if(!legend.horiz) ncoll <- 1

	if(pExpl==1) {
		nlevels.i <- length(unique(expl.data[,1]))
		if( !miss.border & nlevels.i > 7 ) {
		# Make sure the factor levels will be visible or give a warning.
		# if(length(pch) != N & length(col)!= N)
		#  warning("Factor levels will not be displayed automatically, as there are too much levels.
		#          Please use the 'col' and 'pch' arguments to make factors visible.")
		}
	} else {
		nlevels.i <- numeric(pExpl)
		
		for(iexpl in 1:pExpl){
			nlevels.i[iexpl] <- length(unique(expl.data[,iexpl]))
		}
	}

	for(iexpl in 1:pExpl) {
		#levelsi  <- levels(expl.data[,iexpl])
		levelsi <- as.character(unique(expl.data[,iexpl]))
		nlevelsi <- nlevels.i[iexpl]

		mvabnami <- c(rbind(mvabnam, matrix("", ncol=n.vars, nrow=(nlevelsi-1)  )))

		factor.shift <- unclass(expl.data[,iexpl])

		move <- rep(0,times=N)

		uniqu <- unique(factor.shift)

		for(i in 1:length(uniqu)) {
			move[factor.shift==uniqu[i] ] <- i 
		}
	
		moveat <- move[N:1]
		move <- (move-(length(uniqu)+1)/2) * 0.6/length(uniqu)
		moveat <- (moveat-(length(uniqu)+1)/2) * 0.6/length(uniqu)
	
		y.axisi <- y.axis - rep(move, times = n.vars)
		y.axisiat <- y.axisat - rep(moveat, times = n.vars)
		# Use 'moveat' for 'at to have the order of the factor starting at the
		# top and going down
	
		if ( miss.border) {
			if(nlevelsi < 10) {
				#bordr <- c("red", "darkgreen", "orange", "plum", "darkred", "darkblue",
				#		"purple","rosybrown", "black", "green", "hotpink", "gold", "brown",
				#			"lightblue","darkgrey")[nlevelsi:1]
				bordr <- c(1:9)[nlevelsi:1]
				border <-  rep(bordr, times = n.vars)
			} else border <- rep(rainbow(nlevelsi)[nlevelsi:1], times = n.vars)
		}  # Use '[nlevelsi:1]' to get the same colors for the factors as in type="p".

		par(mar=c(4,6,4,2) +0.1, fg=fg, oma=c(1,2,1,1))

		if(legend) {
			# Specify space for the legend.
			if(legend.horiz) {ncoll <- ceiling(nlevelsi/(5+5*row))}
			else ncoll <- min(3, ncoll)
	
			if(!is.null(xlim)) {
				xlimi <- c(xlim[1], xlim[2] + (xlim[2]- xlim[1])*ncoll/6)
			} else xlimi <- xlim
			
		} else xlimi <- xlim

                if (miss.col) colt <- if(type=="bx") "white" else "black" 
                else colt <- col

		do.call( "boxplot", c(list(as.vector(mvabund.object)~y.axisi, xlab="",	horizontal=TRUE,ylab="", main=main, names=mvabnami, at=unique(y.axisiat)*10, las=las,cex.axis=caxis, ylim=xlimi, axes=FALSE, xlim= c(min(y.axisi-0.5), max(y.axisi+0.5))*10, col=colt,border=border), dots))

		# Add a label for the y axis.
#       		mtext(ylab,side=2,line=6,col=colab, cex=clab )
       
		# Add a label for the x axis.
       		mtext(xlab,side=1,line=4,col=colab, cex=clab )

		if (scale.lab=="s" ) {
      			# Adjust some axis details.
      			if (minmva==0)  { posx=-0.05  }
          		else  posx=0
               	} else {
           		seque<-c(minmva ,axTicks(1), mx)
           		sequenc<-c("",as.character(axTicks(1)),"")
           		posx <- minmva
      		}
      
		posy <- 0.5 

      		# Specify below axis,left and Top of plot.
		# Draw the final box around the plot (right edge)
		axis(side=2,at=unique(y.axisiat)*10, labels=c("",mvabnami[1:(length(mvabnami)-1)]),las=las, pos=posx, col=fg, outer=TRUE, cex.axis=caxis)
#		axis( side=1,las=las,pos=min(y.axisi-0.5)*10,col=fg,cex.axis=par("cex.axis"))
      	
#		if  ((transformation!="no") & (t.lab=="o" )) {
#			# Add an additional third axis showing transformations.
#			axis(side=1,las=las, pos=min(y.axisi-0.5)*10, col=fg, at=seque, labels=sequenc, cex.axis=par("cex.axis")) 
#		} else {
			axis( side=1, las=las , pos=min(y.axisi-0.5)*10, col=fg, cex.axis=caxis, at=seque,labels=sequenc)
#		}

#		rec.x = seque[length(seque)] + 1/20*seque[length(seque)] 
#		rect(xleft=posx, ybottom= min(y.axisi-0.5)*10 , xright=rec.x, ytop=max(y.axisi+0.5)*10, border=fg)
		
                rect(xleft=posx, ybottom= min(y.axisi-0.5)*10, xright=tick.max+0.5, ytop=max(y.axisi+0.5)*10, border=fg)
                

		#do.call( "boxplot", c(list(as.vector(mvabund.object)~y.axisi, xlab="" ,
		#		horizontal=TRUE,ylab="", main=main, names=mvabnami, at=unique(y.axisiat)*10,
				# The names are the variable lables for the tickmarks
		#		las=las,cex.axis=0.6, ylim=xlimi,
		#		xlim= c(min(y.axisi-0.2), max(y.axisi+0.2))*10, col=col,border=border), dots))

		# Add a label for the y axis.
		#mtext(ylab,side=2,line=4.5,col=colab,cex=par("cex.lab")*par("cex")*clab )

		# Add a label for the x axis.
		#mtext(xlab,side=1,line=4.5,col=colab,cex=par("cex.lab")*par("cex")*clab )

		#if ((transformation!="no") & ( t.lab=="o" )){
			# Add an additional third axis showing transformations.
		#	axis(side=3,at=axis3, labels=ax3lab, cex.axis=0.6, las=las,col=fg)   
		#}
               
                if (miss.pch) pcht <- 1
                else pcht <- unique(pch)

		if(legend) {
		   tmp    <- suppressWarnings(as.numeric(levelsi))
		   natmp  <- is.na(tmp)
		   leg    <- rep("", times = length(tmp))
		   leg[natmp]   <-  substr(levelsi[natmp], 1, 20)
		   leg[!natmp]  <- zapsmall(as.numeric(levelsi[!natmp]),20)
		   legend(legend=leg, x = "topright", horiz=legend.horiz, col=unique(border)[nlevelsi:1], cex=cex, ncol=ncoll, pch= pcht, text.col="black" ) }
#                   else {
#                       legend(title=substr(factornames[iexpl],1,20), legend=leg, x = "topright", horiz=legend.horiz, col=unique(border)[nlevelsi:1], cex=cex, ncol=ncoll, pch= unique(pch),text.col="black" ) }
	# Use '[nlevelsi:1]' to get the right colors in the legend.
		
	}
############  END draw a boxplot  #################

############  Start a scatterplot  #################
} else {
    sh <- .1  # adjust shifting space
    if(!legend.horiz) ncoll <- 1
    if(pExpl==1) {
	nlevels.i <- nlevels(expl.data[,1] )
    } else {
	nlevels.i <- numeric(pExpl)
	for(iexpl in 1:pExpl)
           nlevels.i[iexpl] <- nlevels(expl.data[,iexpl])
    }

#    par(mfrow=c(1, 2))
    par(mgp=c(1.4,1,2), mar= c(3,8,1,1), oma= c(.5,.5,.5,.1) )
    for(iexpl in 1:pExpl) {            
 	#levelsi  <- levels(expl.data[,iexpl])
	levelsi <- as.character(unique(expl.data[,iexpl]))
	nlevelsi <- nlevels.i[iexpl]

	factor.shift <- unclass(expl.data[,iexpl])

	move <- rep(0,times=N)

	uniqu <- unique(factor.shift)

	for(i in 1:length(uniqu)) move[factor.shift==uniqu[i] ] <- i 
	
	move <- (move-(length(uniqu)+1)/2) * 0.6/length(uniqu)
	
	y.axisi <- y.axis - rep(move, times = n.vars)

	if(miss.pch ) {
	    if(nlevelsi <=7) {
		pchr <- c(1:5, 15,17)[1:nlevelsi]
		pch <- rep(0, times = N)
		for(lev in 1:nlevelsi) pch[factor.shift == lev] <- pchr[lev]
		pch <-  rep(pch, times = n.vars)
	    } else pch = 1	
        } 
	if (miss.col) {
	    if (nlevelsi<= 9) colr <- c(1:9)[1:nlevelsi]
	    else colr <- rainbow(nlevelsi)
	    col <- rep(0, times = N)
	    for(lev in 1:nlevelsi) 
	 	col[factor.shift == lev ] <- colr[lev]
	    col <-  rep(col, times = n.vars)
	}
		
	########### factor plot #################
	# Shift overlapping points.
	if (shift)   {
	    y.axisi <- y.axisi+shiftpoints(y.axisi,c(mvabund.object), sh=sh, method=2)
	  # message("Overlapping points were shifted along the y-axis to make them visible.")
	}

        xlimi <- xlim	   
#        par(mgp=c(1.4,1,2), mar= c(3,6,1,.5),oma= c(.5,.5,.5,.1) )
	if(legend==TRUE) {
	    # Specify the space for the legend.
	    if(legend.horiz) { 
		ncoll <- ceiling(nlevelsi/(5+5*row))
	    } else ncoll <- min(3, ncoll)
	
	    if(!is.null(xlim)) 
		xlimi <- c(xlim[1], xlim[2] + (xlim[2]- xlim[1])*ncoll/6)
	} 

	plot(mvabund.object,y.axisi,ylab="",xlab="", main="", ylim=c(0,n.vars+0.5), pch=pch, type=type, xlim=xlimi, col=col, cex=cex, lwd=lwd, axes=FALSE)

       # add title 	
       mtext(main, side=3, line=-1.2, cex=clab, font=2, adj=0.5)
#       mtext("(b)", side = 3, cex = 1.8, at=-1.8, line=-1.2)

	# Specify some axis details.
	if (scale.lab=="s" ){
	    if (minmva==0)  posx=-0.05 
            else  posx=0
	} else {
	    posx=minmva
	}

#        seque<-c(minmva ,axTicks(1), mx)       
        seque <- axTicks(1)
        if (missing(xaxis.labels)) 
            sequenc<-as.character(axTicks(1))        
        else sequenc <- xaxis.labels

	posy <- 0.5
        # 'at' is turned around to adjust for that the yaxis data is starting
	# at top of y-axis going downwards,
	axis( side=2, at=c((n.vars+0.5), n.vars:1,0.5),	labels=c("" , substring(mvabund.colnames,first=1,last= max(10, nchar(mvabund.colnames, type="char") )  ),""), las=las,pos=posx, col=fg, outer=TRUE, cex.axis=caxis)
        # A label for the y axis.
#	mtext(ylab,side=2,line=6,cex=clab,col=colab)
	# A label for the x axis.
	mtext(xlab,side=1,line=2,cex=clab,col=colab)
        axis(side=1,las=las,pos=posy,col=fg,cex.axis=caxis,at=seque,labels=sequenc)

#	rec.x = seque[length(seque)] + 1/20*seque[length(seque)] 
#	rect(xleft=posx, ybottom= posy , xright=rec.x, ytop=(n.vars+0.5), border=fg)
	rect(xleft=posx, ybottom= posy, ytop=(n.vars+0.5), xright=tick.max+0.5, border=fg)

	if(legend) {
#	   if (is.null(leg)) {
                if (is.character(levelsi)) {
                    leg <-  substr(levelsi, 1, 12)
		}    
                else {
                    tmp <- suppressWarnings(as.numeric( levelsi ))
                    natmp <- is.na(tmp)
                    leg[natmp] <-  substr(levelsi[natmp], 1, 10)
                    leg[!natmp] <- zapsmall(as.numeric(levelsi[!natmp]),10)
                }
#	   }
	   par(xpd=NA)
	   legpch <- unique(pch)
           legcol <- unique(col)	   
           legend("bottomright", inset=c(.15,0.08), legend=leg, horiz=legend.horiz,col=legcol, cex=0.9*cex, ncol=ncoll, pch=legpch, bg=0, lwd=lwd, lty=rep(0,length(leg)), x.intersp=0.01)
	}
    } # for(iexpl in 1:pExpl) 
} ## end scatterplot - if (type="bx")

if(n.vars < p) {
    if(miss.varsubset) 
        tmp <- " \n(the variables with highest total abundance)"   
    else {tmp <- " (user selected)" }

    tmp2 <- paste(colnames(mvabund.object), collapse = ", ")
    message("Only the variables ", tmp2, " were included in the plot", tmp, ".")
}

if(!any(is.na(subset))) 
    message("Only the subset ", allargs$subset,	"  of the cases was included in the plot(s) (user selected).")

}

