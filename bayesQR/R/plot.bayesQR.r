plot.bayesQR <- function(x, var=NULL, quantile=NULL, burnin=0, credint=c(.025,.975), plottype=NULL,  
                         main=NULL, xlab=NULL, ylab=NULL, xlim=NULL, ylim=NULL,...){

  # Error handling
  pandterm <- function(message) {
      stop(message, call. = FALSE)
  }

  if (length(plottype)!=1){
	  pandterm("Plottype should be 'quantile' or 'trace' or 'hist'")
	} else if (!(plottype %in% c("trace","quantile","hist"))){
	  pandterm("Plottype should be 'quantile' or 'trace' or 'hist'")
	}

  # Number of estimated quantiles
	nqr <- length(x)

	# Number of variables in each quantile regression
  nvar <- length(x[[1]]$names)

	# If no variable is specified, plot all variables 
	if (is.null(var)){
	  var <- 1:nvar

	# If variable name is given, find index
	} else if (is.character(var)){
	  if (!all(var %in% x[[1]]$names)){
		  pandterm("Variable name does not exist in x")
		}
	  var <- which(x[[1]]$names %in% var)

	# If variable index is given, check if exists
	} else if (is.numeric(var)){
	  if (!all(var %in% c(1:nvar))){
		  pandterm("Incorrect variable index")
		}
	}


  # Quantile plot
	#oooooooooooooooooooooooooo
	if (plottype=="quantile") {
	  
    if (nqr<3) {
		  pandterm("To few estimated quantiles to create a quantile plot")
		}	

    # Summarize the x
    QRsumobj <- summary(object=x, burnin=burnin, credint=credint)

    # Create quantile plot(s)
		z1 <- FALSE; z2 <- FALSE
		for (i in 1:length(var)){

      # Create plotdata
      plotdata <- matrix(sapply(QRsumobj,"[[","betadraw"),nrow=nvar)[var[i],]
			plotdata <- cbind(sapply(QRsumobj,"[[","quantile"),matrix(plotdata,ncol=3,byrow=TRUE))

			# Plot limits
			if (is.null(xlim)) xlim <- c(0,1)
			if (is.null(ylim)) ylim <- c(min(plotdata[,2:4]),max(plotdata[,2:4])); z1 <- TRUE

			# Plot labels
			if (is.null(main)) main <- ""
			if (is.null(xlab)) xlab <- "quantile"
			if (is.null(ylab)) ylab <- paste("Beta ",var[i],sep=""); z2 <- TRUE

  		# Plot axes/box in correct scale
  		plot(x=NULL, y=NULL, xlim=xlim, ylim=ylim, main=main, xlab=xlab, ylab=ylab, ...)

      # calculate a small value (outside/below the box region)
      small <- min(plotdata[,2:4])-(max(plotdata[,2:4])-min(plotdata[,2:4]))
  
      # plot credible interval
  		polygon(x=c(plotdata[1,1],plotdata[,1],plotdata[nqr,1]), 
  		        y=c(small,plotdata[,4],small),col="grey",border=FALSE)
  
  		polygon(x=c(plotdata[1,1],plotdata[,1],plotdata[nqr,1]), 
  		        y=c(small,plotdata[,3],small),col="white",border=FALSE)
  
      # plot Bayes estimate
  		points(x=plotdata[,1],y=plotdata[,2],typ="o",lty=2)
  
  		# some aesthetic additions
  		points(x=plotdata[,1],y=plotdata[,3],typ="l",col="darkgrey")
  		points(x=plotdata[,1],y=plotdata[,4],typ="l",col="darkgrey")
  		box(lwd=1.3,col="white")
  		box(lwd=1.3,col="black")

			# ask for user input to go to next plot
			if (i < length(var)){
        ans <- readline("Do you want to see the next plot (type 'y' or 'n'):\n")
        while ((ans != "y") & (ans != "n")) ans <- readline("Incorrect input, type 'y' or 'n':\n")
        if (ans == "n") break 
			}

			# set some changed parameters back to input values
			if (z1) ylim <- NULL
			if (z2) ylab <- NULL
		}
	}

  # Trace plot
	#oooooooooooooooooooooooooo
	if (plottype=="trace") {

	  # if quantiles are specified, check if x
		if (!is.null(quantile)){
  		allquant <- sapply(x,"[[","quantile")
  		if(!all(quantile %in% allquant)){
 		    pandterm("Specified quantile does not exist in x")
		  }
		  loopvec <- which(allquant %in% quantile)
		} else {
		  loopvec <- 1:nqr
		}

    # set 'ans'
		ans <- "n"

    # Loop trough quantiles
		for (i in loopvec){
		
		  # Loop trough all specified variables
			for (ii in var){
				plotdata <- x[[i]]$betadraw[,ii]
        plot(plotdata[(burnin+1):length(plotdata)], typ="l", xlab="iteration", ylab="beta", 
				     main=paste("Quantile: ", x[[i]]$quantile, " - Beta ", ii))

        # Ask user input
  			if (!((i==tail(loopvec,n=1))&(ii==tail(var,n=1)))){
          ans <- readline("Do you want to see the next plot (type 'y' or 'n'):\n")
          while ((ans != "y") & (ans != "n")) ans <- readline("Incorrect input, type 'y' or 'n':\n")
          if (ans == "n") break 
  			}
			}
    if (ans == "n") break 
		}
	}

  # Histogram plot 
	#oooooooooooooooooooooooooo
	if (plottype=="hist") {

	  # if quantiles are specified, check if x
		if (!is.null(quantile)){
  		allquant <- sapply(x,"[[","quantile")
  		if(!all(quantile %in% allquant)){
 		    pandterm("Specified quantile does not exist in x")
		  }
		  loopvec <- which(allquant %in% quantile)
		} else {
		  loopvec <- 1:nqr
		}

    # set 'ans'
		ans <- "n"

    # Loop trough quantiles
		for (i in loopvec){
		
		  # Loop trough all specified variables
			for (ii in var){
				plotdata <- x[[i]]$betadraw[,ii]
        hist(plotdata[(burnin+1):length(plotdata)], breaks=100, prob=TRUE, xlab="beta", 
				     main=paste("Quantile: ", x[[i]]$quantile, " - Beta ", ii))

        # Ask user input
  			if (!((i==tail(loopvec,n=1))&(ii==tail(var,n=1)))){
          ans <- readline("Do you want to see the next plot (type 'y' or 'n'):\n")
          while ((ans != "y") & (ans != "n")) ans <- readline("Incorrect input, type 'y' or 'n':\n")
          if (ans == "n") break 
  			}
			}
    if (ans == "n") break 
		}
	}
} 
