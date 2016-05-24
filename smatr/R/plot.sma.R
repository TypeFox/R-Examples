
plot.sma <- function(x, which=c("default","residual","qq"),  use.null=FALSE, add=FALSE, type='o', 
	xaxis=NULL, yaxis=NULL, xlab=NULL, ylab=NULL, pch=NULL, col=NULL, lty=NULL, from=NULL, to = NULL, log=x$log, 
	frame.plot = TRUE, tck=par("tck"),p.lines.transparent=NA, ...){

	# function used to make colours transparent alpha = 0 means fully transparaent
	make.transparent <- function(col, alpha=1) {
  		tmp <- col2rgb(col)/255
	rgb(tmp[1,], tmp[2,], tmp[3,], alpha=alpha)
	}

	#preprocessing ------------------------------------------------------------
	obj <- x  # this needed for consistency with plot
	if(obj$gt == "none"){
		ngrps <- 1	
	}
	else{
		groups <- levels(obj$data[,3])
		ngrps <- length(groups)
	}
	
	whichplot <- match.arg(which)
	
	#---colors--------------------------------
	#user-defined colors
	if(!is.null(col[1])){	
		if(length(col)== 1 &&  ngrps > 1)
			col<-rep(col[1],ngrps); #check right vector length 
	} else {
	#default colors
		col <- c("blue2",  "goldenrod1", "firebrick2", "chartreuse4", "deepskyblue1", "darkorange1", 
		"darkorchid3", "darkgrey", "mediumpurple1", "orangered2", "chocolate", "burlywood3",
		"goldenrod4", "darkolivegreen2", "palevioletred3", "darkseagreen3", "sandybrown", "tan", 
		"gold", "violetred4", "darkgreen")
		col <- rep(col, ceiling(ngrps/length(col)))
	}
	
	#---symbols--------------------------------	
	#user-defined symbols
	if(!is.null(pch[1])){	
		if(length(pch) == 1 && ngrps > 1)
			pch <- rep(pch[1],ngrps) #check right vector length 
	} else #default SYMBOLS		
		pch <- rep(1,ngrps) 
	#---line type--------------------------------	#user-defined symbol	
	if(!is.null(lty[1])){
			if(length(lty) == 1 && ngrps > 1)
			lty<-rep(lty[1],ngrps) #check right vector length 
	} else #default SYMBOLS
		lty <- rep("solid", ngrps)
	#-----------------------------------------------------------------------------
	# DATA PLOT	
	if(whichplot == "default"){

		 #obtain data--------------------------------
		 Y <- obj$data[,1] 
		 X <- obj$data[,2] 
	    
	    #log trasnformations--------------------------------
	    log_data <- obj$log	#logging of data based on transformation applied in sma. Allows scaling of axes ot be changed, while maintaing correct transformation of data
		XLog <- YLog <- 0
		if((log_data == "x") | (log_data == "xy")){ XLog=1; X = 10^X}
		if((log_data == "y") | (log_data == "xy")){ YLog=1; Y = 10^Y}

		#axis labels--------------------------------
	    #determine axis labels if not already specified
	    if(is.null(xlab)){
	    	xlab <- names(obj$data)[2]
			if(XLog)
				xlab <- paste(xlab, "[log scale]")
		}
	    if(is.null(ylab)){
	    	ylab <- names(obj$data)[1]
			if(YLog)
				ylab <- paste(ylab, "[log scale]")
		}

		#SETUP AXES--------------------------------
		if(add==FALSE)
		{
			#use nice plotting if appropriate
			if(!is.null(xaxis)  && !is.null(yaxis)){
				
				#Deteremine axis limits if missing - caluclated on transformed data. 				#add 5% white space around plots, then back transform if necessary
				if (is.null(xaxis$limits)){
					Range_x <-range(obj$data[,2])
					Buffer <- (Range_x[2]-Range_x[1])*0.05
					xaxis$limits <- c(Range_x[1]-Buffer, Range_x[2]+Buffer) 
					if(XLog)xaxis$limits <- 10^xaxis$limits
				}
				if (is.null(yaxis$limits)){
					Range_y <-range(obj$data[,1])
					#add 4% white space around plots (like R default in plot)
					Buffer <- (Range_y[2]-Range_y[1])*0.04
					yaxis$limits <- c(Range_y[1]-Buffer, Range_y[2]+Buffer)   
					if(YLog) yaxis$limits <- 10^yaxis$limits
				}
				
				#Make plot
				nicePlot(xaxis,yaxis,log=log,xlab=xlab, ylab=ylab, 
					frame.plot = frame.plot, tck=tck,...)
			}
			else
				plot(X,Y, type='n', log=log, xlab=xlab, ylab=ylab,
					frame.plot = frame.plot, tck=tck,...)
		}
	
		#add datapoints	--------------------------------
		if(type %in% c("o","b", "p")){
			if(obj$gt == "none")
				points(X, Y, col = col[1], pch=pch[1],...)
			else{
				for(i in 1:ngrps){
					iref  <- as.character(obj$data[,3]) == groups[i]
					points(X[iref], Y[iref], col =col[i], pch=pch[i],...)
				}
			}
		}
		
		#add lines --------------------------------
		if(type %in% c("o","b", "l")){
			
			#decide which coefficients to use: alternative (default) or null
			if(use.null==FALSE)
				coef <- obj$coef
			else
				coef <- obj$nullcoef
			
			#determine end points for lines
			if(is.null(from[1])){  #based on fitted values
				for(i in 1:ngrps){
					from[i] <- as.numeric(obj$from[i])  
					to[i] <- as.numeric(obj$to[i])
				}
			} else {  #user defined
				if(length(from) == 1){
					from <- rep(from[1], ngrps)
					to <- rep(to[1], ngrps)
				}
			}
			
			#add lines to plot
			for(i in 1:ngrps){
				#coefficients
				a <- coef[[i]][1,1]
				B <-  coef[[i]][2,1]
				
			    p <- obj$groupsummary$p[i]

				if(!is.na(p.lines.transparent))
					col.tr <- make.transparent(col[i], max(0, (1 - p/p.lines.transparent)))
				else
					col.tr <-  col[i]

				#choose line according to log-trsnaformation used in fitting data, even if different transformation used for axes
				if(log_data=="xy")
	        		curve(10^a*x^B, from[i], to[i], add=TRUE,col = col.tr, lty= lty[i],...)
    	   	 	if(log_data=="x")
        			curve(a+B*log10(x), from[i], to[i], add=TRUE, col = col.tr, lty= lty[i],...)
        		if(log_data=="y")
        			curve(exp(a+x*B), from[i], to[i], add=TRUE, col = col.tr, lty= lty[i],...)
       	    	if(log_data=="")
        			curve(a + x*B, from[i], to[i], add=TRUE,  col = col.tr, lty= lty[i],...)
        	}
        }
	}

	# RESIDUAL PLOT	
	if(whichplot == "residual")
	{
		 #obtain data--------------------------------
		Y <- fitted.sma(obj, type = "residuals")
		X <- fitted.sma(obj, type = "fitted")

		
		#axis labels--------------------------------
	    #determine axis labels if not already specified
	    if(is.null(xlab)) xlab <- paste("Fitted values (",names(obj$data)[2], " v ",
			names(obj$data)[1],")")  
	    if(is.null(ylab)) ylab <- paste("Residual values (",names(obj$data)[2], " v ",
			names(obj$data)[1],")")  

		#SETUP AXES--------------------------------
		if(!add){ #use default plotting options 
			plot(X,Y, type='n', xlab=xlab, ylab=ylab, frame.plot = frame.plot,...)
   		}
		
		#add datapoints	--------------------------------
		if(type %in% c("o","b", "p")){
			if(obj$gt == "none")
				points(X, Y, col = col[1], pch=pch[1],...)
			else{
				for(i in 1:ngrps){
					iref <- as.character(obj$data[,3]) == groups[i]
					points(X[iref], Y[iref], col =col[i], pch=pch[i],...)
				}
			}
		}
	}

	# QQ PLOT	
	if(whichplot == "qq")
	{
		 #obtain data--------------------------------
		Y <- fitted.sma(obj, type = "residuals")
		
		#axis labels--------------------------------
	    #determine axis labels if not already specified
	    if(is.null(xlab)) xlab <- "Normal quantiles"
	    if(is.null(ylab)) ylab <- paste("Residual values (",names(obj$data)[2], " v ",
			names(obj$data)[1],")")  

		#SETUP AXES--------------------------------
		if(add==FALSE){ #use default plotting options 
			qqnorm(Y, xlab=xlab, ylab=ylab,...)
			qqline(Y)
   		}
	}

}
	
	
