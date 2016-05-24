################################################################################
# PLOT.MVABUND: Plot functions for mvabund objects (multivariate abundance data)
################################################################################
default.plot.mvabund <- function(x, 
				y, 
				type="p", 
				main="", 
				xlab="", 
				ylab="",
  				col= if(type=="bx") "white" else "black", 
				fg= "grey", 
				pch=1, 
				las=1,
  				write.plot="show", 
				filename="plot.mvabund", 
				n.vars= 12,
  				overall.main, 
				var.subset=NA, 
				subset=NA, 
				transformation="log", 
				scale.lab="ss",
   				t.lab="o", 
				mfrow=if(!two.objects | (type=="bx")) { 1 } 
				else if(write.plot=="show") min(25, n.vars) 
				else min(9,n.vars), 
				mfcol=NULL,
  				shift=TRUE, 
				border=	if(two.objects & length(col)==1) { c("red","blue")}else "black", 
				checks = TRUE, 
				add.line=FALSE, 
				line.col="black",
  				keep.window=FALSE, 
				ask=TRUE, ...) 
{ 

 dev <- dev.list()
 dev.name <- getOption("device")
 #if (!is.null(dev)) dev.off() # close previous window

 if(is.null(dev.name))
 	stop("Make sure that the 'device' option has a valid value, e.g. 'options(device = 'windows')'.	Allowed values here are 'windows', 'win.graph', 'x11', 'X11'.") 

 ########## BEGIN check transformation, t.lab, scale.lab ###########

 if(substr(transformation, 1,1)=="n" | transformation==""){
    transformation <- "no"
 } else if(substr(transformation, 1,1)=="l"){ transformation <- "log"
 } else if(substr(transformation, 1,1)=="s" & transformation =="sqrt"){
    transformation <- "sqrt"
 } else if(substr(transformation, 1,1)=="s" & transformation =="sqrt4"){
    transformation <- "sqrt4"
 } else stop("You have passed an invalid 'transformation'")
 # scale.lab <- substr(scale.lab,1,1)
 # if(!scale.lab %in% c("r", "s")) stop("You have passed an invalid 'scale.lab'")
 t.lab <- substr(t.lab,1,1)
 if(!t.lab %in% c("o", "t")) stop("You have passed an invalid 't.lab'")

# allargs = all arguments that are passed to the function
allargs <- match.call(expand.dots=FALSE)  # argument list
dots <- allargs$...
if(!is.null(dots$log)){
    dots$log <- NULL
    warning("argument 'log' not implemented in 'plot.mvabund'")
}

allargs[[1]] <- NULL
allargs$x <- NULL
allargs$y <- NULL
allargs$... <- NULL
args <- lapply(allargs, eval, parent.frame()) # value of arguments
targs <- match.call(call = sys.call(which = 1), expand.dots = FALSE)
# targs = plot(x = tasm.cop ~ treatment, col = as.numeric(block)) the original input

 ########## BEGIN mvabund data check, pipe to alternate Function call if needed 
if (missing(x)) { stop("The mvabund object 'x' is missing.") }

 if (missing(y)) two.objects <- FALSE
 else {   # Find out which object is mvabund and which function needs to be called
    two.objects <- TRUE

    if (checks) { 
       mvabund.colnames <- colnames(as.matrix(unabund(x)))
       mvabund.colnames2 <- colnames(as.matrix(unabund(y)))
       if(any(c(is.null(mvabund.colnames), is.null(mvabund.colnames2)))) ch <- TRUE
       else ch <- sapply(1:length(mvabund.colnames), 
                  function(x) mvabund.colnames[x] == mvabund.colnames2[x])
       if(!all(ch)) stop("The two mvabund objects 'x' and 'y' seem not to consist of the same variables in the same order.")
    } 

    if (is.data.frame(x) & !is.factor(y) ) {
         # x is not mvabund and y is a data.frame.
         # as they are both data frames, but the response in a formula cannot
         # be a data frame, there is no interpretation as formula possible.
         stop("The mvabund object 'x' is missing.")
    }     
    else if (is.mvabund(y) & is.mvabund(x)) {  
        allargs$n.vars <- min(n.vars, NCOL(x))
        mvabund.object.2 <- as.matrix(unabund(y))
        if (any(!is.na(list(subset)))) 
           mvabund.object.2 <- mvabund.object.2[c(subset),]          
    } 
    else if(is.mvabund(x) & is.factor(y)) {	
        # x is mvabund data and y factor design
        allargs$n.vars <- min(n.vars, NCOL(x))
        expl.data <- as.data.frame(y)
         # Check if all columns in x are factors.
        fac <- rep(TRUE, (ncol(expl.data)))
        for (i in 1:ncol(expl.data)) 
             if (!is.factor(expl.data[,i])) fac[i] <- FALSE                
        if (any(!fac)) warning("Only the factor variables ", paste((1:ncol(expl.data))[fac], collapse =", "), " of x will be plotted.")
        if(shift) message("Overlapping points were shifted along the y-axis to make them visible.")
        cat("\n PIPING TO 1st MVFACTOR \n")
        do.call("default.plotMvaFactor", c(list(x=x, y=expl.data), allargs, dots))
        return(invisible())
    }
    else if(is.mvabund(y) & is.factor(x)) {
        # y is assumed to be mvabund and x assumend to be factors.
        allargs$n.vars <- min(n.vars, NCOL(y))
        expl.data <- as.data.frame(x)
        # Check if all columns in x are factors.
        fac <- rep(TRUE, (ncol(expl.data)))
        for(i in 1:ncol(expl.data)) 
            if (!is.factor(expl.data[,i])) fac[i] <- FALSE                    
        if(any(!fac)) warning("Only the factor variables ", paste((1:ncol(expl.data))[fac], collapse =", "), " of x will be plotted.")
        if(shift) message("Overlapping points were shifted along the y-axis to make them visible.")
        cat("\n PIPING TO 2nd MVFACTOR \n")
        do.call("default.plotMvaFactor", c(list(x=y, y=expl.data), allargs, dots))
        return(invisible())
     }
     else if (!is.mvabund(x) & !is.data.frame(y)){ 
        foo <- eval(targs$x)
        if( inherits(foo, "formula")) {
            cat("\n PIPING TO 1st MVFORMULA \n")	
            do.call("default.plot.mvformula", foo, allargs, dots)
            return(invisible())
        }  
        else {
            # If x is not a factor, y is assumed to be the response and the
            # mvabund object, pass on to plot.mvformula,
            # Define the formula.mva <- y ~ x.
            formula.mva <- as.formula(paste(deparse(targs$y,width.cutoff = 500),
          	"~", deparse(targs$x,width.cutoff = 200), sep=""))

            if(transformation != "") { 
               if (min(x,na.rm=TRUE )< 0 | min(y,na.rm=TRUE) < 0)
                   stop("no negative values allowed in the data if a transformation is used")
               if (transformation == "log") { 
                   miny <- eval(min(y[y!=0],na.rm=TRUE))
                   minx <- eval(min(x[x!=0],na.rm=TRUE))
                   formula.mva <- update( formula.mva, log(. +miny)/miny ~ log(. + minx)/minx)
               } else if(transformation == "sqrt4") { 
                   formula.mva <- update( formula.mva, (.)^0.25 ~ sqrt(sqrt(.)) )	
               } else if(transformation == "sqrt") { 
                   formula.mva <- update( formula.mva, sqrt(.) ~ sqrt(.) )
               }
            }
            
            allargs$n.vars <- min(NCOL(y), 12)
             
            cat("\n PIPING TO 2nd MVFORMULA \n")            	
	    do.call("default.plot.mvformula", c(list(x=formula.mva), allargs, dots))
            return(invisible())
        }
    } 
   else if (!is.mvabund(y) & !is.data.frame(x)){
      	# x is assumed to be mvabund.
        foo <-  eval(targs$y, envir =  parent.frame())
     	if( inherits(foo, "formula") ) {
            data <- eval( targs$x, envir =  parent.frame())
            plot.mvformula(x=foo, y=data, allargs, dots)
            return(invisible())
        }
       	# x is the mvabund and the response, pass on to plot.mvformula.
        else {
      	# Define the formula.mva <- x ~ y
            formula.mva <- as.formula(paste(deparse(targs$x,width.cutoff = 200),
      	        	"~", deparse(targs$y,width.cutoff = 500), sep=""))
            if(transformation != "") {
       	       if( min(x,na.rm=TRUE )< 0 | min(y,na.rm=TRUE) < 0)
                  stop("no negative values allowed in the data if a transformation is used")
               if(transformation == "log") { 
                  minx   <- eval(min(x[x!=0],na.rm=TRUE))
                  miny   <- eval(min(y[y!=0],na.rm=TRUE))
                  formula.mva <- update(formula.mva, log(. +minx)/minx ~ log(. + miny)/miny )
	       } else if(transformation == "sqrt4") { 
                  formula.mva <- update( formula.mva, (.)^0.25 ~ sqrt(sqrt(.)) )
	       } else if(transformation == "sqrt") { 
                  formula.mva <- update( formula.mva, sqrt(.) ~ sqrt(.) )
               }
            }

            allargs$n.vars <- min(n.vars, NCOL(x))
            cat("\n PIPING TO 3rd MVFORMULA \n")
       	    do.call( "default.plot.mvformula", quote=FALSE, c(list(x=formula.mva),allargs,dots))
            return(invisible())
    	} # if (class(foo)==formula)
    } #  (!is.mvabund(y) & !is.data.frame(x))
    else
       stop("unknown data type of x and y")
  } # if (!missing(y) | !is.null(y))
 #### END mvabund data check, FACTOR & FORMULA PLOTS should have been piped elsewhere!! ##

 ########## PREPERATIONS FOR PLOT, OUTPUT, EXTRA INPUT, N, P, SUBSET ######### 
 if (write.plot!="show")   {	
	if (write.plot=="eps" | write.plot=="postscript") {
      		postscript(paste(filename,".eps", sep=""))
    	} else if (write.plot=="pdf") {
		pdf(paste(filename,".pdf", sep=""))
    	} else if (write.plot=="jpeg" ) {
		jpeg(paste(filename,".jpeg", sep=""))
    	} else if (write.plot=="bmp" ){
		bmp(paste(filename,".bmp", sep=""))
    	} else if (write.plot=="png" ){
		png(paste(filename,".png", sep=""))
	}
    	# Specify the window where to draw the plot.
    	dev.curr <- dev.cur()
    	on.exit(dev.off(which = dev.curr))
 }

 miss.varsubset <- missing(var.subset) | is.null(var.subset) | !is.numeric(var.subset)
 # Change logical var.subset to numerical var.subset, if necessary. Note that
 # NA values are logical as well, but should be excluded here.
 if(!miss.varsubset){
 	if(is.logical(var.subset) & any(!is.na(list(var.subset)))) var.subset <- which(var.subset[!is.na(list(var.subset))])
 }
 if(length(dots)>0) { 
 # Delete arguments in ... that are defined lateron and cannot be used twice
 # in the plot function.
 	if(type=="bx") { 
               deactive <- c("xlim", "ylim", "axes", "horizontal", "cex.axis")
    	} else deactive <- c("xlim", "ylim", "axes", "horizontal", "names", "cex.axis") 

        deactivate <- (1:length(dots))[names(dots) %in% deactive ]
        for (i in length(deactivate):1) {
        	dots[ deactivate[i] ]<-NULL }				#fixed up deactivate, caused compile error (v2.10)

        	dots <- lapply(dots, eval, parent.frame(), parent.frame())
        	
		if ("cex.lab" %in% names(dots)) clab <- dots$cex.lab
        	else clab <- 1
        
		if ("col.lab" %in% names(dots)) colab <- dots$col.lab
        	else colab <- par("col.lab")
 } else {
 	clab <- 1
        colab <- par("col.lab")
 }

 mvabund.object.1 <- as.matrix(unabund(x))

 # Initiate n.vars before deleting x.
 if(!any(is.na(subset))) mvabund.object.1 <- mvabund.object.1[c(subset),, drop=FALSE]
  
 N <- nrow(mvabund.object.1)     # number of sites
 p <- ncol(mvabund.object.1)     # number of organism types
 
 if (any(c(p,N)==0)) stop("The mvabund object 'x' has invalid dimensions.")
 
 var.subset<-as.vector(var.subset) 

 ########## PREPERATIONS FOR PLOT, OUTPUT, EXTRA INPUT, N, P, SUBSET #########
  
 ############################################################################
 ## option a: one mvabund object is passed,                                 #
 ## a scatterplot of abundance vs. variable number is produced              #
 ############################################################################
 
 if( !two.objects ) {
        n.vars <- min(n.vars, p)
        cat("Kicking off BoxPlot sequence \n")
 	rm(x)
    	mvabund.colnames <- colnames(mvabund.object.1)

    	if(is.null(mvabund.colnames)) mvabund.colnames <- paste("Variable", 1:p)

    	opp <- par("ask", "col.lab", "cex.lab", "mar", "fg", "mgp" ,"mfcol" ,"mfrow")

    	if(keep.window & write.plot=="show") opp$mfrow <- opp$mfcol <- NULL
    	if(!is.null(mfcol)) mfrow <- mfcol else opp$mfrow <- opp$mfcol <- NULL
    
	# Make sure that not a new window is drawn with the next plot after
    	# exiting the function.
    	if (length(mfrow)==1)  {
        	columns <-ceiling(sqrt(mfrow))
            	row     <-columns-1
            	if (row*columns<mfrow) row<-columns
            	mfrow = c(row,columns)
         	if(write.plot=="show") {
         		if(all(opp$mfrow==c(row,columns))) opp$mfrow <- opp$mfcol <- NULL
         	}
        
	} else {
         	row <- mfrow[1]
         	columns <- mfrow[2]
        }

    	if(write.plot=="show") on.exit(par(opp))
    	# Upon exiting the function, reset all graphical parameters to its value
    	# at the beginning.
       
    	if((write.plot=="show" & is.null(dev) ) & (!is.null(mfrow))) {
        	if (mfrow[2] > mfrow[1]){
            		width   <- 8
            		height  <- max(mfrow[1]*width/mfrow[2],5)*1
        	} else {
            		height  <- 8
            		width   <- max(height*mfrow[2]/mfrow[1],4)*1
        	}
     	
		# Close the old window before proceeding.
        	dev.off()
 
        	do.call(dev.name, args=list(height=height,width=width))
        	dev.flag <- TRUE
    
	} else dev.flag <- FALSE

       	if (!is.null(mfcol)) par(mfcol=mfrow)    
       	else if (!is.null(mfrow)) par(mfrow=mfrow) 
   
         ######### BEGIN edit var.subset, n.vars and mvabund.objects  #########
        # subset allows double variables
        var.subset.dim <- length(var.subset)
        if (!is.numeric(var.subset)) {
             if (n.vars > p)
        	stop("You have passed an invalid number of variables 'n.vars' to be included in the plot.")

             sum.mvabund.object.1 <- t(mvabund.object.1)%*%matrix(1,ncol=1,nrow=N)
        	
	     # Find abundance ranks OF MVABUND.OBJECT.1.
             var.subset <- order(sum.mvabund.object.1, decreasing = TRUE)
        
	     # Ensure no more than n.vars in var.subset.
             if (n.vars < length(var.subset)) var.subset <- var.subset[1:n.vars]
        	var.subset.dim<-length(var.subset)

        # Arrange data to plot requested var.subset (default - n.vars most abund).
      	} else if (p<max(var.subset))  {
      	     stop("You have passed an invalid 'var.subset'")
        # Do some dimension checks for the subset.
      	} else if (n.vars!=var.subset.dim) { 
	     n.vars<- var.subset.dim
	}
      
  	mvabund.object.1 <- mvabund.object.1[,var.subset, drop=FALSE]
  	mvabund.colnames <- mvabund.colnames[var.subset]
  
  	######### END edit var.subset, n.vars and mvabund.objects #########
   	# Get variable numbers, the abundances are plotted against them
   	# yaxis data should start at the top of y-axis going downwards.
   	y.axis <- rep(n.vars:1 , each=N)

  	# Get max value for y axis before any transformations.
   	mx <- max(mvabund.object.1, na.rm=TRUE)
   	scale.lab <- substr(scale.lab, start=1, stop=1)

	#Get Tickmarks using function and transformation type!
	xmin <- min(mvabund.object.1[mvabund.object.1>0])
	xmax <- max(mvabund.object.1)
	xticks <- axisTicks(transformation, xmax, xmin, t.lab)

  	if (scale.lab=="s") {
  	    # Do a 's'tandard plot: 0 is included in the x axis
     	    #xlim <- c(-0.1,max(mx+mx/10, 0.1 ))
     	    # Don't allow scale.lab="s" when there are negative values in the data
     	    # (as log transf is used).
     	    if ( min(mvabund.object.1,na.rm=TRUE) < 0 )
        	stop("'scale.lab' cannot be 's' if the data contains negative values")
  	} else if ( scale.lab=="r") {
  		# Do an 'r' plot: R's default is used
      		#xlim <- NULL
  	} else stop ("undefined value for 'scale.lab'")


  	if( transformation!="no" ) {
	     if (t.lab=="o" ) {
    		if (scale.lab=="s") ylim <- c(1, max(mx+mx/10, 2 ) ) else ylim <- NULL

    		if(dev.flag) do.call(dev.name, args=list(height=height,width=width)) 
		else do.call(dev.name, args=list())
        
    		if(any(c(mvabund.object.1)!=0) ) {
        	     suppressWarnings( plot( c(mvabund.object.1),type="n",axes=FALSE, 
			xlab="", ylab="",ylim=NULL,log="y"))
        		#axis3 <- axTicks(2)
    		} else {
                     plot( mvabund.object.1,y.axis,type="n",axes=FALSE,xlab="",
          		ylab="",xlim=xlim,ylim=c(0,n.vars+0.5))
       			# Get the transformation-labels for third axis.
       			#axis3 <- axTicks(2)
       			#lengthax3 <- length(axis3)
       			#if (lengthax3 > 5) axis3 <- axis3[-(5:(lengthax3-1))]
    		}
    			
		dev.off()
    		#ax3lab <-  as.character(axis3) 
    		#if (scale.lab=="s") {
    		#ax3lab <- c("0", ax3lab )
    		#axis3 <- c(0,axis3) }
    	     }

	     ######### BEGIN transformation #########
  		# Transform data, if required.
  	     if (transformation=="log") {
     	         if (max(mvabund.object.1,na.rm=TRUE)==0)
        	     stop("The mvabund object 'x' only consists of zero-abundances")
		 minNon0 <- min(mvabund.object.1[mvabund.object.1!=0],na.rm=TRUE)  
    		 mvabund.object.1 <- log(mvabund.object.1+minNon0)-log(minNon0)    
    			
		 # plot title: 'Abundances [log(y/min+1) scale]'
    		 transf.lab <- expression(paste("Abundances  ", 
		 bgroup("[", paste(log, bgroup( "(",frac(y,min)+1,")"), " scale"),"]")))
		 #if (t.lab=="o" ) {axis3 <- log(axis3+minNon0)-log(minNon0)}
  	     } else if (transformation=="sqrt4") {
     		 mvabund.object.1 <- (mvabund.object.1)^0.25
     		 # plot title:  'Abundances [y^{0.25} scale]'
     		 transf.lab <- expression(paste("Abundances  ",
        			bgroup("[", paste(sqrt(y,4)," scale") ,"]")))
     			
		 #if (t.lab=="o") {axis3<-(axis3)^0.25 } 
   		
	     } else if (transformation=="sqrt") {	
		  mvabund.object.1 <- sqrt(mvabund.object.1)	
		  # plot title: 'Abundances [y^{0.5} scale]'
     		  transf.lab <- expression(paste("Abundances  ",
        			bgroup("[", paste(sqrt(y)," scale") ,"]")))
		  #if (t.lab=="o") { axis3 <- sqrt(axis3)  }
    	     }
	     ######### END transformation #########
    	} else transf.lab <- 'Abundances'
   
 	if( write.plot!= "show")  dev.set(which = dev.curr)
    	# Specify the window where to draw the plot. 	
 
 	######### BEGIN some calculations for better axis scaling #########
 	if (missing(main)) main <- "Abundances" 
 	if (missing(xlab)) xlab <- transf.lab
 	if (missing(ylab)) ylab <- "Species"
   
 	# Get min for the correction of posx.
 	minmva <- min(mvabund.object.1,na.rm=TRUE)
	
	# Get max value for y axis after transformation.
 	#mx <- max(mvabund.object.1,na.rm=TRUE)

   
 	#if (scale.lab=="s" ) {
 	
		#xlim <- c(-0.1,mx+mx/10)
    	
		#if (mx == 0) {  # i.e. all data ==0
     			#mx <- 0.1
     			#xlim <- c(-0.1, 0.1)
   	
		#} else if (mx < 0.8) {
     			#potl<- -floor(log(mx)/log(10) )
     			#mx <- ceiling(mx *(10^potl)) /(10^potl)
   		#} else  if (mx>6){
    	
			#if (mx<10) potl <- - ceiling(log(mx)/log(10) ) 
			#else potl <- - floor(log(mx)/log(10) )
      			
			#mx <- ceiling(mx *(10^potl)) /(10^potl)
    		#} else mx<-ceiling(mx)
   
		#seque <- seq(from=0, to=mx, by=mx/10)
		#reset xlim, to final limit!   		
		#xlim <- c(-0.1,mx+mx/10)

		#if (mx > 1000) {
      			#potence <- 10^(floor(log(mx)/log(10)))
      			#sequenc <- as.character( round(seque/potence,digits=2 ))
      			#xlab <- paste( xlab, " in ", potence)
   		#} else if (transformation!="no") {
        		#sequenc<- as.character(round(seque,digits = 2 ) )
	  	#} else sequenc <- as.character(round(seque,digits = 2 ) )

 	#} else {seque <- NULL
      		#sequenc <- NULL
      	 	#labels <- NULL
      	 	#xlim <- NULL }

	#Create a Robust xlim, based on Tick Generation
	xlim <- c(-0.1,1.05*max(mvabund.object.1,na.rm=TRUE))

 	######### END some calculations for better axis scaling #########

 	######### BEGIN plot #########
 	# Draw either a BOXPLOT
 	if (type=="bx") {
 	   if(!is.null(dots$boxwex)) { 
		boxwex.i <- dots$boxwex
          	dots$boxwex <- NULL
      	   } else boxwex.i <- 0.7
      
	   if(!is.null(dots$names)) {  
		names.i <- dots$names
      		dots$names <- NULL
      	   } else {
        	names.i <- mvabund.colnames[n.vars:1]
        		# n.mn <- nchar(names.i, type="chars")
        		# if(any(n.mn > 12))  names.i[n.mn > 12] <-
        		# paste(substr(names.i[n.mn > 12], 1, 7),"...\n...",
        		# substr(names.i[n.mn > 12], n.mn[n.mn > 12]-7, n.mn[n.mn >12]) )
                names.i <- substring(names.i,first=1,last= max(10,nchar(names.i, type="chars")))
           }

	   # Create Outer Margin for Optimal Layout
       	   par(mar=c(4, 6, 3, 1)+0.1, fg=fg, oma=c(1,2,1,1)) 
       	   do.call( "boxplot", c(list(as.vector(mvabund.object.1)~y.axis, xlab="" ,
      	  	horizontal=TRUE, ylab="", main=main , names=names.i,
       	 # 'n.vars:1' because data is starting at top of y-axis going downwards
       	# names are the variable lables for the tickmarks
       		las=las,cex.axis=0.8, ylim=xlim,col=col,border=border, boxwex=boxwex.i, axes=FALSE), dots))
	   # Add a label for the y axis.
       	   mtext(ylab,side=2,line=5,col=colab, cex=par("cex.lab")*par("cex")*clab )
	   # Add a label for the x axis.
       	   mtext(xlab,side=1,line=2.5,col=colab, cex=par("cex.lab")*par("cex")*clab )

	   if (scale.lab=="s" ) {
      		# Adjust some axis details.
      		if (minmva==0)  { posx = -0.05  }
          	else  posx=0
           } else {
           	#seque<-c(minmva ,axTicks(1), mx)
           	#sequenc<-c("",as.character(axTicks(1)),"")
           	posx <- minmva
      	   }
	   posy <- 0.5 

      	   # Specify below axis,left and Top of plot.
	   # Draw the final box around the plot (right edge)
	   rect(xleft=posx, ybottom= posy , xright=xlim[2], ytop=(n.vars+0.5), border=fg)
	   axis(side=2,at=c((n.vars+0.5), n.vars:1,0.5), labels=c("",names.i,""), las=las, pos=posx, col=fg, outer=TRUE, cex.axis=0.75)
       	
	   # Draw additional third axis showing transformations.
	   if ((transformation!="no") & ( t.lab=="o" )){
       		axis(side=1,pos=posy, at=xticks$x.tic, labels=xticks$x.ticlab, cex.axis=0.75, las=las,col=fg)
           } else {
		axis(side=1, las=las , pos=posy , col=fg, cex.axis=0.75, at=xticks$x.tic, labels=xticks$x.ticlab)
	   }  
       
 	} else { # type=="bx") 
	# Draw a DOTPLOT.
          	sh <- 0.5
    
    		# Shift overlapping points.
    		if (shift) {
         		y.axis <- y.axis+shiftpoints(y.axis,c(mvabund.object.1), sh=sh, method=2)
       			message("Overlapping points were shifted along the y-axis to make them visible.")
    		}
      
    		par(mgp=c(1.5,1,0),mar=c(4, 6, 4, 2) + 0.1)
    	
		#This is the main plot call!
        	cat("\n\nABOUT TO PLOT THE FUNCTION \n\n")

		do.call( "plot", c(list(mvabund.object.1,y.axis,ylab="",xlab="", main=main,ylim=c(0,n.vars+0.5), pch=pch, axes=FALSE,  type=type, xlim=xlim,col=col),dots))

      		if (scale.lab=="s" ) {
      		    # Adjust some axis details.
      		     if (minmva==0)  { posx=-0.05  }
          	     else  posx=0
               	} else {
           	     #seque<-c(minmva ,axTicks(1), mx)
           	     #sequenc<-c("",as.character(axTicks(1)),"")
           	     posx <- minmva
      		}
      
		posy <- 0.5 
#browser()
      		# Specify the left axis.
		if (is.null(dots$yaxt)) {
	      	      axis( side=2,at=c((n.vars+0.5), n.vars:1,0.5), 
		      labels=c("" , substring(mvabund.colnames,first=1,
     	              last= max(10, nchar(mvabund.colnames, type="char"))),""),
        	      las=las,pos=posx, col=fg, outer=TRUE, cex.axis=0.6)
                } 
		# Add a label for the y axis.
      		mtext(ylab,side=2,line=4,cex=par("cex.lab")*par("cex")*clab, col=colab )
      
		# Add a label for the x axis.
      		mtext(xlab,side=1,line=2.5,cex=par("cex.lab")*par("cex")*clab, col=colab )
      
		# Draw additional third axis showing transformations.
		if ((transformation!="no") & ( t.lab=="o" )){
       		    axis(side=1,pos=posy, at=xticks$x.tic, labels=xticks$x.ticlab, cex.axis=0.75, las=las,col=fg)
        	} else {
		    axis(side=1, las=las , pos=posy , col=fg, cex.axis=0.75, at=xticks$x.tic, labels=xticks$x.ticlab)
		}
      
		# Draw the final box around the plot (right edge)
		rect(xleft=posx, ybottom= posy , xright=xlim[2], ytop=(n.vars+0.5), border=fg)
 	} # type=="bx") 

	if(n.vars < p) {
 	    if(miss.varsubset) 
                 tmp <- " \n(the variables with highest total abundance)"   
	    else { 
                 tmp <- " (user selected)" }
                 message("Only the variables ", paste(colnames(mvabund.object.1), collapse = ", ")," were included in the plot", tmp,".")
 	}
   
 	if(!any(is.na(subset))) message("Only the subset ", allargs$subset, "  of the cases was included in the plot(s) (user selected).")

 	######### END plot #########
  } else { # !two.objects

        cat("\n PIPING to TWO object Plots \n")
        # Draw plot for two.objects.
  	if  ((N!=NROW(mvabund.object.2))| (p!=NCOL(mvabund.object.2)))
        	stop("the dimesions of 'x' and 'y' do not match")

  	#############################################################################
  	## If 2 mvabund objects with the same dimensions are passed                 #
  	## plots of one against the other are produced,                             #
  	## with the n.vars (default 12) most abundant variables plotted.            #
  	#############################################################################
    	if (checks)  {
            mvabund.colnames <- colnames(mvabund.object.1)
            mvabund.colnames2 <- colnames(mvabund.object.2)
            if (any(c(is.null(mvabund.colnames), is.null(mvabund.colnames2)))) ch <- TRUE 
	    else ch <- sapply(1:length(mvabund.colnames), function(x) mvabund.colnames[x] == mvabund.colnames2[x])
        	
	    if(!all(ch))
           	stop("The two mvabund objects 'x' and 'y' seem not 
		to consist of the same variables in the same order.")
    	 }
    
    	# Make labels for x and y axis.
    	# Use the first element, because if the string is longer
    	# than 200 it is assumed to be a vector of character values
    	xlabel <- deparse(substitute(x),  width.cutoff = 200)[1]
    	ylabel <- deparse(substitute(y),  width.cutoff = 200)[1]

    	# When plot.mvabund is called from boxplot the name is the whole object,
    	# change the name in that case.
    	if(substr(xlabel, 1,9) == "structure")    xlabel <- "x"
    	if(substr(ylabel, 1,9) == "structure")    ylabel <- "y"
    	rm(x)
    	rm(y)

    	# Get title for boxplot.
    	mainbx <- "Abundances"
    	if (! is.null(colnames(mvabund.object.1)) ) names <- colnames(mvabund.object.1) 
        else if (! is.null(colnames(mvabund.object.2)) ) names <- colnames(mvabund.object.2) 
        else names <- paste ("Species", 1:p)    

    	if (missing(main)) { 
        	main <- paste( "\n",substring(names,first=1, last= max(8, max(nchar(names, type="char")))) )
    	} else if (length(main)==1) { 
        	mainbx <- main
        	main <- rep(main, times=p) 
    	} else if (is.numeric(var.subset) & length(main)==length(var.subset)) {
        	msav <- main
        	main <- rep.int(NA, times=p)
        	main[var.subset] <- msav
        } else if (length(main)!=p)
      		stop("the length of 'main' does not match to the dimension of response variables")

    	dataAll <- rbind(mvabund.object.1,mvabund.object.2)

    	######### BEGIN edit var.subset, n.vars and mvabund.objects #########
        var.subset.dim <- length(var.subset)
        if (!is.numeric(var.subset)) {
        
        	if (n.vars > p) stop("You have passed an invalid number of 
						response vectors to be included in the plot.")
  		# Find abundance ranks.
		sum.dataAll <- t(dataAll) %*% matrix(1,ncol=1,nrow=2*N)
                var.subset <- order(sum.dataAll, decreasing = TRUE)
           
		# Ensure no more than n.vars in var.subset.
           	if (n.vars < length(var.subset)) var.subset <- var.subset[1:n.vars]
            	
		var.subset.dim <- length(var.subset)
           	# Arrange data to plot requested var.subset (default - n.vars most abund).
        } else if (p < max(var.subset)) {
                stop ("You have given an invalid 'var.subset'")
        
	} else if (n.vars!=length(var.subset)) {
        	# Do some dimension checks for subset & n.vars.
        	n.vars<- var.subset.dim
        }   
        
        mvabund.object.1 <- mvabund.object.1[,var.subset, drop = FALSE] 
        mvabund.object.2 <- mvabund.object.2[,var.subset, drop = FALSE]
        dataAll <- rbind(mvabund.object.1,mvabund.object.2)
        main    <- main[var.subset]
        names   <- names[var.subset]
    	######### END edit var.subset, n.vars and mvabund.objects #########

    	######### BEGIN establish row and column sizes #########   
    	if (all(is.null(c(mfcol, mfrow))))  {
        	perwindow <- par("mfrow")
        	# Open new window if none is open yet.
                mfr <- TRUE
    	} else if (is.null(mfcol)){
        	perwindow <- mfrow
        	mfr <- TRUE
    	} else {
        	perwindow <- mfcol
        	mfr <- FALSE
    	}

    	if (length(perwindow)==1) {
        	windows<-ceiling(n.vars/perwindow)  # number of windows
        	perwind<-min(n.vars,perwindow)      # number of plots per window
        
		if (prod(par("mfrow"))>=perwind) {
            		row     <- par("mfrow")[1]
            		columns <- par("mfrow")[2]
        	} else {
            		columns <-ceiling(sqrt(perwind))
            		row     <-columns-1
            		if (row*columns<perwind) row<-columns
        	}
    	
	} else {
        	row     <- perwindow[1]
          	columns <- perwindow[2]
          	perwindow <- row*columns
          	windows<-ceiling(n.vars/perwindow)
          	perwind <- min(n.vars,perwindow)
    	}

    	opp <- par("ask", "col.lab", "cex.lab", "mar", "fg", "mgp" ,"mfcol" ,"mfrow")
    	if((keep.window | all(opp$mfrow==c(row,columns))) & write.plot=="show")
    	# Avoid that a new window is drawn with the next plot after exiting the function.
    	opp$mfrow <- opp$mfcol <- NULL
    
    	# Reset all graphical parameters to its value at the beginning, upon
    	# exiting the function.
    	if(write.plot=="show") on.exit(par(opp))
    
    	if(write.plot=="show" & is.null(dev) ) {
        	if (columns > row){
            		width   <- 10
            		height  <- max(row*width/columns,5)
        	} else {
            		height  <- 10
            		width   <- max(height*columns/row,4)
        	}
        
		# Close the old window before proceeding.
        	dev.off()
	        do.call(dev.name, args=list(height=height,width=width))
        	dev.flag <- TRUE

    	} else dev.flag <- FALSE


	# Reduce the Margins of the ith Plot & Add an Outer Margin for main title
	par(mar=c(3,3,3,1),oma=c(1,1,2,1))

	######### END establish row and column sizes #########

    	xlim <- NULL
    	mxAll <- ceiling(max(dataAll,na.rm=TRUE))

    	scale.lab1 <- substr(scale.lab, start=1, stop=1)
    	if ( scale.lab1=="s"  ) {
      		mxi <- numeric(n.vars)
     		for (i in 1:n.vars) {
        		# Get max value for y axis before transformation.
        		mxi[i] <- ceiling(max(dataAll[,i],na.rm=TRUE))
      		}
      		wh.max <- which( abs( mxi - max(mxi, na.rm=TRUE)) < 1e-06  )[1]
    	}
	
    	if (transformation!="no") {
    		# Calculate transformations and transformed tickmarks + ticklabels.
          	if (type=="bx" & t.lab=="o") {

      			# Get max value for y axis before transformation.
      			# minim <- min(dataAll,0, na.rm=TRUE)
      			if ( scale.lab1=="s" ) xlim=c(1, max(mxAll+mxAll/30,1))
      
      			if(dev.cur() == 1) {
         			if(dev.flag) do.call(dev.name, args=list(height=height,width=width)) 
				else do.call(dev.name, args=list())
      			}
	
      			if (any(c(dataAll)!=0) ) {
        			suppressWarnings(plot(c(dataAll),type="n",axes=FALSE, 
        					xlab="", ylab="",xlim=xlim, ylim=xlim, log="y")) 
		
			} else  plot(c(dataAll),type="n",axes=FALSE, 
            				xlab="", ylab="",xlim=xlim, ylim=xlim) 

      			axis3 <- axTicks(2)
      			ax3lab <- as.character(axis3)
      		
			if (scale.lab1=="s" | min(dataAll,na.rm=TRUE)==0) {
      				# Negative values are not allowed if a transformation is used.
          			axis3 <- c(0, axis3)
          			ax3lab <- c("0", ax3lab) 
          		}
          
			axis4 <- ax4lab <- NULL
       			dev.off()    # close the window
      
		} else if (t.lab=="o" ) {
       			# Calculate ideal axis before transformation.
      			axis3 <- ax3lab <- axis4 <- ax4lab <- list()

      			if(dev.flag) do.call(dev.name, args=list(height=height,width=width)) 
			else do.call(dev.name, args=list())

      			# for "ss" use the same for all variables
      			if(scale.lab=="ss") loopind <- wh.max 
			else loopind <- 1:n.vars

	      		for (i in loopind) {
      				# Get max value for y axis before transformation.
      				# mxi <- ceiling(max(dataAll[,i],na.rm=TRUE))
      				# minim <- min(dataAll[,i],0, na.rm=TRUE)
      				if ( scale.lab1=="s"  ) xlim=c(1, max(mxi[i]+mxi[i]/30,1))
      
      				plot(c(mvabund.object.1[,i], mvabund.object.2[,i]),type="n",axes=FALSE, 
         				xlab="", ylab="",xlim=xlim,ylim=xlim)
      			
				maxmva1 <- maxmva2 <- NULL    
        
      				if(any(c(mvabund.object.1[,i], mvabund.object.2[,i])!=0 )) {
      					maxmva1 <- axTicks(1)[length(axTicks(1))]
      					maxmva2 <- axTicks(2)[length(axTicks(2))]
      					suppressWarnings(plot(c(mvabund.object.1[,i], mvabund.object.2[,i]),
        							type="n",axes=FALSE,xlab="", ylab="",
									xlim=xlim,ylim=xlim, log="xy"))
        			} 
      
				#define a sequence to extact only half of the tix marks
				xtick.seq <- seq(1,length(axTicks(1)),by=2)
				ytick.seq <- seq(1,length(axTicks(2)),by=2)
				
				axis3[[i]]<-c(axTicks(1)[xtick.seq], maxmva1)
      				ax3lab[[i]]<- as.character(axis3[[i]])
      				axis4[[i]]<- c(axTicks(2)[ytick.seq],maxmva2) 
      				ax4lab[[i]]<- as.character(axis4[[i]])
      
      				if (scale.lab1=="s" | min(mvabund.object.1[,i],na.rm=TRUE)==0) {
	      				# negative values are not allowed if a transformation is used
          				axis3[[i]] <- c(0, axis3[[i]])
          				ax3lab[[i]] <- c("0", ax3lab[[i]]) 
          			}
      
				if (scale.lab1=="s" | min(mvabund.object.2[,i],na.rm=TRUE)==0) {
      					# negative values are not allowed if a transformation is used
                			axis4[[i]] <- c(0, axis4[[i]])
          				ax4lab[[i]] <- c("0", ax4lab[[i]])
        			} 
      			}

      			if(scale.lab=="ss") {
        			for(i in 1:n.vars) {
          				axis3[[i]]  <- axis3[[wh.max]]
          				ax3lab[[i]] <- ax3lab[[wh.max]]
          				axis4[[i]]  <- axis4[[wh.max]]
          				ax4lab[[i]] <- ax4lab[[wh.max]]
        			}
      			}
      
      			dev.off() # close the window 
      		}

    		############ BEGIN transform data, if required and get axisticks  ##########
    		############ and axislabels for transformed axis                  ##########

		#Get appropriate axes
		tick.min <- min(dataAll[dataAll!=0],na.rm=TRUE)
		tick.max <- max(dataAll,na.rm=TRUE)
		tick <- axisTicks(transform=transformation, max=tick.max, min=tick.min, tran.lab=t.lab)


    		if (transformation=="log")  {  
cat("Performing Log Transformation \n")    
        		minNon0 <- min(dataAll[dataAll!=0],na.rm=TRUE)    
        		mvabund.object.1<-log(mvabund.object.1+minNon0)-log(minNon0)   
        		mvabund.object.2<-log(mvabund.object.2+minNon0)-log(minNon0)
        		transf.lab <-  transfy.lab <- '[log(y+1) scale]'
	      		transf.xlab <- '[log(x+1) scale]'
      
			dataAll <- rbind(mvabund.object.1,mvabund.object.2)
      	
			if (t.lab=="o")  {
      				if(type=="bx") axis3 <- log((axis3) + minNon0)-log(minNon0) 
				else {
       					for (i in 1:n.vars){
       						# Calculate transformation-axis if chosen with t.lab
             					axis3[[i]]<-log((axis3[[i]])+minNon0)-log(minNon0)
       						axis4[[i]]<-log((axis4[[i]])+minNon0)-log(minNon0)
					}
				}
       				# Adjust cex: smaller for transformation-axis
       				cex.axis<-0.7
     			} else cex.axis <- 0.9

     		} else if (transformation=="sqrt4") {
cat("Performing 4th Root Transformation \n")     
       			mvabund.object.1<-mvabund.object.1^0.25
       			mvabund.object.2<-mvabund.object.2^0.25
       			transf.lab  <- expression(paste(bgroup("[",paste(sqrt(y,4)," scale") ,"]")))
       			transf.xlab <-  "[x^0.25 scale]"
	     		transfy.lab <-  "[y^0.25 scale]"
	     		dataAll <- rbind(mvabund.object.1,mvabund.object.2)
       
			if (t.lab=="o")  {
       				if(type=="bx") { axis3<-(axis3)^0.25 } 
				else {
       					for (i in 1:n.vars) {
       						axis3[[i]]<-(axis3[[i]])^0.25
       						axis4[[i]]<-(axis4[[i]])^0.25  
					}
				}
       				cex.axis<-0.7
       			} else cex.axis<-0.9

     		} else if (transformation=="sqrt") {
cat("Performing SQRT Transformation \n")     
       			mvabund.object.1<-sqrt(mvabund.object.1)
       			mvabund.object.2<-sqrt(mvabund.object.2)
       			transf.lab  <-  expression(paste(bgroup("[",paste(sqrt(y)," scale") ,"]")))
	     		transf.xlab <- "[sqrt(x) scale]"
	     		transfy.lab <- "[sqrt(y) scale]"
       			dataAll<-rbind(mvabund.object.1,mvabund.object.2)
       
			if (t.lab=="o" ) {
       				if(type=="bx") { axis3<-sqrt(axis3) } 
				else {
       					for (i in 1:n.vars) {
       						axis3[[i]]<-sqrt(axis3[[i]])
       						axis4[[i]]<-sqrt(axis4[[i]]) 
					} 
				}
       				cex.axis<-0.7 
			} else cex.axis <- 0.9
    		}
     
    		if(type != "bx" & t.lab=="o") {
     			for (i in 1:n.vars){
    				# If the last value is very close to the second last, it can be skipped,
    				# Alternative: skip the second last value
        			laxis3 <- length(axis3[[i]])
        			laxis4 <- length(axis4[[i]])
        			last3need <- ((axis3[[i]])[laxis3]-(axis3[[i]])[laxis3-1]) <
          					((axis3[[i]])[laxis3]-(axis3[[i]])[1])/10
        			
				last4need <- ((axis4[[i]])[laxis4]-(axis4[[i]])[laxis4-1]) <
          					((axis4[[i]])[laxis4]-(axis4[[i]])[1])/10

        			if(last3need) {
            				axis3[[i]] <- (axis3[[i]])[-laxis3]
            				ax3lab[[i]]<- (ax3lab[[i]])[-laxis3]
        			}
        
				if( last4need) {
            				axis4[[i]] <- (axis4[[i]])[-laxis4]
            				ax4lab[[i]]<- (ax4lab[[i]])[-laxis4]
        			}
     			}
    		}
 	} else { transf.lab  <-  transf.xlab <- transfy.lab <- "Abundances"
      		 cex.axis <-  0.9 }

     	if (missing(overall.main)) overall.main <- transf.lab

    	if (missing(xlab))  { 
        	if(type=="bx") xlab <- transf.xlab
        	xlab <- paste(substr(xlabel, 1, 8), transf.xlab)
    	}
    
	if (missing(ylab))  { 
    		if (type=="bx")  ylab <- "Species" 
		else ylab <- paste(substr(ylabel,1,8), transfy.lab)
    	}
    	######### END transformation #########
   	
	if( write.plot!= "show")  {
    		# Draw the plot to the specified device.
    		dev.set(which = dev.curr)
    	}
    
    	######### BEGIN plot #########
    	if (write.plot=="show" & n.vars>perwindow & type!="bx" & ask) par(ask=TRUE)
    	if ( scale.lab1=="r"  ) { xlim <- NULL }

    	###### BEGIN Window LOOP #########
        if(all(is.null(c(mfcol, mfrow))) & windows==1) {
    		# Leave the window as it is, but if a new one needs to be drawn
    		# (not enough space on window as planned), ask
        	par(ask=TRUE)
    	} else if (mfr) { 
		par(mfrow=c(row,columns))
    		# Might produce an error if perwindow was chosen too large.
       	} else { 
		par(mfcol=c(row,columns))
    	}
     
	if(type == "bx")  {
cat("Beginning Two Object Boxplot\n")
		dev.off()
		stop("\nERROR: A boxplot is not applicable for paired data, please use \nanother plot type\n")
#        		if ( t.lab=="o" & scale.lab1=="s") {
#        			# Get max value for y axis before transformation.
#             		mxAll <- ceiling(max(dataAll,na.rm=TRUE))
#             		xlim=c(0,mxAll+mxAll/30)
#         	}
#        		if(!is.null(border)) border <- border[length(border):1]
#        		if(!is.null(col)) col    <- col[length(col):1]
#        
# 		par(mar=c(6, 5.5, 5.5, 2) +0.1, fg=fg) 
#        		namesbx <- rep("",times = 2* n.vars)
#        		namesbx[ (2*(1:n.vars))] <- paste("\n",
#           	substring(names,first=1, last= min(12, nchar(names, type="char"))) )
#        
# 		# Get variable numbers, the abundances are plotted against them.
#        		y.axis <- rep(2*(1:n.vars) , each=N)
#        		y.axis <- c(y.axis+0.4, y.axis-0.4)
#        
#        		if(!is.null(dots$main)) { 
# 			main.i <- dots$main
#        			dots$main <- NULL
#        		} else main.i <- mainbx
#        
# 		if(!is.null(dots$names)) { 
# 			names.i <- dots$names
#        			dots$names <- NULL
#        		} else {
#        			names.i <- namesbx  
# 		}
# 
#        		if(!is.null(dots$at)) { 
# 			at.i <- dots$at
#        			dots$at <- NULL
#        		} else at.i <- c(- 0.4 + rep(2*(n.vars:1), each=2)+rep(c(0, 0.8),times=n.vars) )
# 
#        		if(!is.null(dots$cex.axis)) { 
# 			cex.axis.i <- dots$cex.axis
#        			dots$cex.axis <- NULL
#        		} else cex.axis.i <- 0.6
#         
#        		do.call( "boxplot", c(list(c(as.vector(mvabund.object.1),
#           			as.vector(mvabund.object.2))~y.axis,
#           				xlab="" , horizontal=TRUE, ylab="", main=main.i , names = names.i,
#        					# names are the variable lables for the tickmarks
#         					las=las,cex.axis=cex.axis.i, xlim=c(1,2*n.vars+1), ylim=xlim,col=col,
#         						border=border, at =at.i), dots))
#         
#        		leg <- rep(par("usr")[1],times=2* n.vars)
#        	
# 		if (length(border)==2) colpoints <- border 
# 		else if (length(col)==2) {
#           		colpoints <- col
#        		} else  colpoints <- c("red", "blue")[c(2,1)]
#        
#        		text.col <- colpoints
#        		points(leg, c(2*(1:n.vars) -0.4,2*(1:n.vars) +0.4 ),
#           			col=rep(c(colpoints),each=n.vars),pch=3 )
#           
#        		# Add a label for the y axis.
#        		mtext( ylab,side=2,line=4.5,col=colab, cex=par("cex.lab")*par("cex")*clab )
#        
#        		# Add a label for the x axis.
#        		mtext( xlab,side=1,line=4.5,col=colab, cex=par("cex.lab")*par("cex")*clab )
#        
# 		# Additional third axis showing transformations.
# 		if ((transformation!="no") & ( t.lab=="o" )){
#        	        	axis(side=3,at=axis3, labels=ax3lab, cex.axis=0.6, las=las,col=fg)
#        		}
#        
# 		# Add a Legend to the plot
# 		legend( x="bottomleft", legend=c(xlabel, ylabel), col=colpoints[c(2,1)],
#           			pch=3, horiz=TRUE, text.col=text.col[c(2,1)], cex=0.8,bty="n" )

 	} else {
     		# Don't allow scale.lab1="s" when there are negative values in the data
     		# (as log transf is used).
     		if ( scale.lab1=="s" & min(dataAll,na.rm=TRUE) < 0 )
        		stop("'scale.lab' cannot be 's' if the data contains negative values")

     		for (l in 1:windows){

     			kchoice <-((l-1)*perwind+1):((l-1)*perwind+perwind)
     			kchoice <- kchoice[kchoice<(n.vars+1)]

     			# The following might produce an error if perwindow was chosen too large.
    			if (mfr & l>1) par(mfrow=c(row,columns))
         		else if (l>1) par(mfcol=c(row,columns))
    			# it's possible to still plot on the first window after the end of the function
    			# only if l >1 --> if only one window
    
    			###### BEGIN Plot LOOP #########
    			mxAll <- ceiling(max(dataAll,na.rm=TRUE))

    			for (i in kchoice) {
    
        			if(scale.lab == "ss") {
           				mxdata <- mxAll
        			} else {
          				# mndata	<-min(dataAll[,i],na.rm=TRUE)
          				mxdata <- max(dataAll[,i],na.rm=TRUE)
        			}
    
				######### BEGIN some scaling calculations #########
      				if ( scale.lab1=="s"  ) {
      
         				if(mxdata == 0 ){
            					mx <- 1
         				} else if (mxdata<0.8) {
         				# Calculate max of plot scale.lab for standard scale.lab.
              					potl 	<- -floor(log(mxdata)/log(10) ) 
              					mx 	<- ceiling(mxdata *(10^potl)) /(10^potl)
          				} else if (mxdata>6){
              					if (mxdata<10) potl	<- - ceiling(log(mxdata)/log(10) ) 
						else potl	<- - floor(log(mxdata)/log(10) )
              					mx 	<- ceiling(mxdata *(10^potl)) /(10^potl)
          				} else mx<- mxAll #mx <- ceiling( max(dataAll[,i],na.rm=TRUE) )

         				seque	<- seq(from=0, to=mx, by=mx/5)  # for scale.lab1 = "s"
         				sequenc <- as.character(seque)

         				# Use standard labels.
         				if ((transformation=="no") | ( t.lab=="t" )) {
             					if (mx > 1000) {
             						potence <- 10^(floor(log(mx)/log(10)))
             						sequenc <- as.character( round(seque/potence,digits=2 ))
             						ylab 	 <- paste(ylab, " in ", potence )
             					} else if (transformation!="no") {
             						sequenc <- as.character(round(seque,digits = 2 ) ) 
						} else	sequenc <- as.character(round(seque,digits = 2 ) ) 
					}

         				sequey 	<- seque
        	 			sequency<- sequenc
					xlim <-c(0, 31 *max(mxdata)/30)
      				}   

    				# Do scaling calculations for R scaling.
    				if (( t.lab=="o") & (transformation!="no"))  {
       
          				seque 	<- axis3[[i]]
          				sequenc	<- ax3lab[[i]]
          				sequey	<- axis4[[i]]
          				sequency<- ax4lab[[i]]
          				
					if ( scale.lab1 =="s"  )
            					xlim <-c(0, 31 *max(seque)/30)
       				}
       
				if (perwindow==1 ) main[i]=paste("\n\n",main[i])
      
				# if (shift) {   # don't use shift in scatterplots!
       				# mvabund.object.2[,i] <- mvabund.object.2[,i] +
       				#   shiftpoints(mvabund.object.2[,i], mvabund.object.1[,i])
       				#}
				xlim <-c(0, 31 *max(tick$x.tic)/30)

      				do.call( "plot", c(list(mvabund.object.1[,i], mvabund.object.2[,i],
        				main= main[i], xlab="",ylab="", xlim=xlim, ylim=xlim ,pch=pch,
        				axes=FALSE ,  type=type,col=col, frame.plot=TRUE, fg=fg), dots))

        			if(add.line) {
            				line.stop <- max(mvabund.object.1,mvabund.object.2,na.rm=TRUE)
            				#lines(c(0,line.stop),c(0,line.stop),col=line.col)
					lines(c(0,xlim),c(0,xlim),col=line.col)
        			}
        
				if ( scale.lab1=="r" & ((t.lab=="t") | (transformation=="no"))  ) {
              				seque <- axTicks(1)
        				if (min(mvabund.object.1[,i],na.rm=TRUE)<=0 & all(0!=seque)) {
            					seque <- c(min(mvabund.object.1[,i],0,na.rm=TRUE), seque  )  
					}
        				sequenc <- substr(as.character(seque), start=1,
          							stop=max(nchar(as.character(axTicks(1)))))
        				sequey	<- axTicks(2)
        
					if (min(mvabund.object.2[,i],na.rm=TRUE)<=0 & all(0!=sequey)) {
            					sequey <- c(min(mvabund.object.2[,i],0,na.rm=TRUE), sequey)  
					}
        				sequency <- substr(as.character(sequey), start=1,
          							stop=max(nchar(as.character(axTicks(2)))))
      				}  

      				###### END scaling calculations ######
				#Use Stephen Calculations
				seque <- tick$x.tic
				sequenc <- tick$x.ticlab

            			# Specify below axis.
      				axis(side=1,at=seque, labels=sequenc ,  las=las, col=fg,  cex.axis=1, pos=NULL)
      
				# Specify left axis.
      				axis(side=2,at=seque, labels=sequenc ,las=las, col=fg, cex.axis=1, pos=NULL)
 
				# Specify a label for the y axis.
				#use modular arithmatic to get label
				ith.yplot <- 0:(n.vars-1)
				is.ylab <- ith.yplot[i]%%columns
				if (is.ylab == 0) ylab.new <- ylab
				else ylab.new <- ""
				mtext(ylab.new,side=2,line=2.5,col=colab, cex=par("cex.lab")*par("cex")*clab*1.3)
      				
				# Specify a label for the x axis.
				#flag if plot is on the last row in window
				ith.xplot <- rep(1:(row*columns),windows)
				is.xlab <- (ith.xplot[i] > row*columns - columns)
				if (is.xlab == TRUE) xlab.new <- xlab
				else xlab.new <- ""
				mtext(xlab.new,side=1,line=2.5,col=colab, cex=par("cex.lab")*par("cex")*clab*1.3)

   			}    
    			###### END Plot LOOP #########

    			if (par("oma")[3] >= 1) mtext(overall.main, outer = TRUE, cex = 1.1* par("cex.main"), col=par("col.main"), line = -0.5)
    		}
 	}
    	###### END Window LOOP #########

    	if(n.vars < p) {
    		if(miss.varsubset) tmp <- " \n(the variables with highest total abundance)"   
		else { tmp <- " (user selected)"  }
      
		message("Only the variables ", paste(colnames(mvabund.object.1),
        			collapse = ", "), " were included in the plot(s)", tmp, ".")
    	}	
    
	if(!any(is.na(subset))) message("Only the subset ", allargs$subset,
        	"  of the cases was included in the plot(s) (user selected).")
    
 }
 ###### END plot  ######
}
###### END FUNCTION CALL ######




