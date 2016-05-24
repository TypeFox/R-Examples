#<<BEGIN>>
hist.mc <- function(x, griddim = NULL, xlab = names(x),ylab = "Frequency", main = "",  ...)
#TITLE Histogram of a Monte Carlo Simulation
#DESCRIPTION
# Shows histogram of a \samp{mcnode} or a \samp{mc} object.
#KEYWORDS hplot
#INPUTS
#{x}<<An \samp{mcnode} or an \samp{mc} object.>>
#[INPUTS]
#{griddim}<<A vector of two integers, indicating the size of the grid of plots. If \samp{NULL}, the grid is calculated to produce a "nice" graph.>>
#{xlab}<<A vector of labels for the x-axis for drawn histograms (those whose \samp{outm(x)!="none"}). May be recycled.>>
#{ylab}<<A vector of labels for the y-axis for drawn histograms. May be recycled.>>
#{main}<<A vector of main title of histograms for drawn histograms. May be recycled.>>
#{\dots}<<Other arguments to be passed to all calls of \samp{hist}.>>
#NOTE
#For Two-dimensional \samp{mc}, the histogram is based on all data (variability and uncertainty) pooled together.
#EXAMPLE
#data(total)
#hist(xVUM3)
#hist(total)

#CREATED 07-08-01
#REVISED 07-08-01
#--------------------------------------------
#
{
# the function beau calculate a nice grid

	beau <- function(n){
		nc <- round(sqrt(n))
		nr <- ceiling(n/nc)
		c(nc,nr)}

    l <- length(x)
    main <- rep(main,l)
	  xlab <- rep(xlab,l)
	  ylab <- rep(ylab,l)


    loutm <- lapply(x,attr,which="outm")
    dimm <- sapply(x,dim)
    n <- sum(dimm[3,]* (loutm=="each") + (loutm!="each" & loutm!="none"))


  	if(is.null(griddim)) griddim <- beau(n)
  	if(prod(griddim) < n) op <- par(mfrow=griddim,ask=TRUE,mar=c(5,4,.2,.2))
    else op <- par(mfrow=griddim,mar=c(5,4,.2,.2))
	
	try({  #to restore par in case of error

    for(i in 1:l){

      if(is.null(loutm[[i]])) loutm[[i]] <- "each"
      if(loutm[[i]][1] == "none") next                                             # Pass outm == none
      for(j in loutm[[i]]){                                                     # loop on the nb of stat, j is the name of the function

        if(j == "each"){
          nvar <- dim(x[[i]])[3]
          if(nvar==1) xlab2 <- xlab[i] else xlab2 <- paste(xlab[i],1:nvar,sep="")
        }
        else {
          func <- get(j,mode="function")                                        # apply the function
          x[[i]] <- apply(x[[i]],c(1,2),func)
          dim(x[[i]]) <- c(dim(x[[i]]),1)
          nvar <- 1                                                             #1 dimension now for this stat
          xlab2 <- paste(j,xlab[i])                                             #change the name with the name of the stat
        }

        if(is.logical(x[[i]])) x[[i]] <- unclass(x[[i]]) * 1                        # unclass to avoid Ops.mcnode

        for(k in 1:nvar){
		      hist(x[[i]][,,k],main=main[i],xlab=xlab2[k],ylab=ylab[i],...)      # loop on nvariates
		    }
      }
    }
  })
  par(op)
	return(invisible())
  }
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#<<BEGIN>>
hist.mcnode <- function(x, ...)
#ISALIAS hist.mc
#--------------------------------------------
{ nom <- deparse(substitute(x))
  x <- list(x)
  names(x) <- nom
  hist.mc(x, ...)}

