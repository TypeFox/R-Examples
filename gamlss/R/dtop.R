# TODO: Add comment
#  It uses the NPL.bands() function of  Ross Darnell from package SMIR
# detrented Owen's plot 
# Author: MS
#########################################################################################
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
#########################################################################################
dtop<-function( object = NULL,
		          xvar = NULL, 
		          resid = NULL,
		    conf.level = c("95", "99"),
		       n.inter = 4,
		   xcut.points = NULL, 
		       overlap = 0,
		    show.given = TRUE, 
		           cex = 1, 
		           pch = 21,
		          line = TRUE,
		           ...
		       )
{
##-------------------------------------------------------------------------------
## this is the NPL.bands() function of  Ross Darnell from package SMIR
##-------------------------------------------------------------------------------	
	NPL.bands <- function (x, conf.level = c("95", "99")) 
	{
		if (!is.numeric(x)) 
			stop("argument must be numeric")
        conf.level <- match.arg(conf.level)
		yn <- table(x)
		yi <- as.numeric(names(yn))
		cn <- as.numeric(cumsum(yn))
		nn <- rep(sum(yn) + 1, length(cn))
		p <- as.numeric(cn/nn)
		if (conf.level == "95") 
		    {
			lambda <- ifelse(nn <= 100, (3.0123 + 0.4835 * log(nn) - 
								0.00957 * (log(nn))^2 - 0.001488 * (log(nn))^3), 
					(3.0806 + 0.4894 * log(nn) - 0.02086 * (log(nn))^2))
		    }
		else 
		    {
			if (conf.level == "99") {
				lambda <- ifelse(nn <= 100, (-4.626 - 0.541 * log(nn) + 
									0.0242 * (log(nn))^2), (-4.71 - 0.512 * log(nn) + 
									0.0219 * (log(nn))^2))
			}
		#else stop("Must be either 0.95 or 0.99")
	    }
		lambda <- sqrt(2 * lambda)
		phi <- pbeta(p, 1/3, 1/3)
		se <- 1/(5.3 * sqrt(nn * (p * (1 - p))^(1/3)))
		phiu <- phi + lambda * se
		phiu <- ifelse(phiu > 1, 1, phiu)
		phil <- phi - lambda * se
		phil <- ifelse(phil < 0, 0, phil)
		pu <- qbeta(phiu, 1/3, 1/3)
		pl <- qbeta(phil, 1/3, 1/3)
		list(x = yi, lower = pl, upper = pu)
	}
##-----------------------------------------------------------------------------
##----------------------------------------------------------------------------- 
	get.intervals <- function (xvar, xcut.points ) 
	{
		if (!is.vector(xcut.points))  {stop(paste("The interval is not a vector."))}
		if ( any((xcut.points < min(xvar)) | any(xcut.points > max(xvar))))
		{stop(paste("The specified `xcut.points' are not within the range of the x: (", min(xvar),
							" , ", max(xvar), ")"))}
		extra <-(max(xvar)-min(xvar))/100000
		int <- c(min(xvar), xcut.points, (max(xvar)+2*extra))
		ii <- 1:(length(int)-1)
		r <- 2:length(int)
		x1 <- int[ii]
		xr <- int[r]-extra
		if (any(x1>xr)) {stop(paste("The interval is are not in a increasing order."))}
		cbind(x1,xr)
	}
##-----------------------------------------------------------------------------
##-----------------------------------------------------------------------------
check.overlap <- function(interval)
{
	if (!is.matrix(interval) ) {stop(paste("The interval specified is not a matrix."))}
	if (dim(interval)[2] !=2) {stop(paste("The interval specified is not a valid matrix.\nThe number of columns should be equal to 2."))}
	crows = dim(interval)[1]
	for (i in 1:(crows-1))
	{
		#if (interval[i,2] != interval[i+1,1]) {interval[i+1,1]=interval[i,2]}
		if (!(abs(interval[i,2]-interval[i+1,1])<0.0001)) {interval[i+1,1]=interval[i,2]}
	}
	return(interval)
}
##-----------------------------------------------------------------------------
##-----------------------------------------------------------------------------
panel.fun<-function (x, y, col = par("col"), bg = NA, pch = par("pch"), 
		cex = par("cex"), col.smooth = "red", span = 2/3, iter = 3, ...) 
{
#	qq <- as.data.frame(qqnorm(y, plot = FALSE))
#	qq$y <- qq$y - qq$x  
#	grid(nx=NA,ny=NA, lwd = 2) 
	xx <- NPL.bands(y, conf.level=conf.level)
	lower <- qnorm(xx$lower)-xx$x      
	upper <- qnorm(xx$upper)-xx$x
	if (!"ylim"%in%names(args))
	{
		points(upper~xx$x, xlab="ordered quantile residuals", 
				ylab= paste(paste(conf.level, "%", sep=""), "confident intevals"),
				bg = "wheat", ylim=c(min(lower[is.finite(lower)]),max(upper[is.finite(upper)])), 
				cex=1,  pch=21,...)
	}
	else
	{
		points(upper~xx$x,xlab="ordered quantile residuals",
				ylab= paste(paste(conf.level, "%", sep=""), "confident intevals"),
				bg = "wheat",  cex=1,  pch=21,...)		
	}    
	points(lower~xx$x, bg = "wheat",  cex=1,  pch=21,...)
	grid(lty = "solid")
	abline(0, 0, lty = 1, col = 2)
	#abline(0, 100000, lty = 2, col = 2)
	if (line)
			{
			fcdf <- ecdf(y) # create a function 
              mm <- fcdf(y) # evaluate it 
          zscore <- qNO(mm) # normalize it 
         dzscore <- (zscore-y) # detrend it
  lines(dzscore[order(y)]~y[order(y)], col=4)
			}

	#no.points <- length(xx$upper)
	#total.points<<-total.points+no.points
}                
## end of local functions here	
##-----------------------------------------------------------------------------
##-----------------------------------------------------------------------------
	# the main function starts here 
	if (is.null(object)&&is.null(resid))  stop(paste("A fitted object with resid() method or the argument resid should be used ", "\n", ""))
	resid <- if (is.null(object)) resid else resid(object)  
	conf.level <- match.arg(conf.level)
	 args <- list(...)
	if(is.null(xvar)) # if there is no x-variable
		{
			   xx <- NPL.bands(resid, conf.level=conf.level)
			lower <- qnorm(xx$lower)-xx$x      
			upper <- qnorm(xx$upper)-xx$x
		if (!"ylim"%in%names(args))
		 {
		   plot(upper~xx$x, xlab="ordered quantile residuals", 
					ylab= paste(paste(conf.level, "%", sep=""), "confident intevals"),
					bg = "wheat", ylim=c(min(lower[is.finite(lower)]),max(upper[is.finite(upper)])), 
					cex=cex,  pch=pch,...)
		 }
		else
		 {
		    plot(upper~xx$x,xlab="ordered quantile residuals",
					ylab= paste(paste(conf.level, "%", sep=""), "confident intevals"),
					bg = "wheat",  cex=1,  pch=21, ...)		
		 }    
			points(lower~xx$x, bg = "wheat",  cex=1,  pch=21,...)
			grid(lty = "solid")
			abline(0, 0, lty = 1, col = 2)
		#	abline(0, 100000, lty = 2, col = 2)
			if (line)
			{
			fcdf <- ecdf(resid) # create a function 
              mm <- fcdf(resid) # evaluate it 
          zscore <- qNO(mm) # normalize it 
         dzscore <- (zscore-resid) # detrend it
  lines(dzscore[order(resid)]~resid[order(resid)], col=4)
			}
	    } 
	else   # if x-variable is defined
	    {
			w  <- if (is.null(object)) rep(1, length(resid)) else  object$weights # geting the weights
			if (all(trunc(w)==w))  xvar <- rep(xvar, w) #
			if(is.null(xcut.points)) 
			{ # getting the intervals automatic
				given.in <- co.intervals(xvar, number=n.inter, overlap=overlap )
				if (overlap==0) given.in <- check.overlap(given.in) 
			}                 
			else
			{ # if xcut.points is set
				given.in <- get.intervals(xvar, xcut.points )
			}
		#	total.points <-0
			coef<-coef1<-NULL
			y <- resid   
			x <- resid
			coplot(y ~ x | xvar, given.values = given.in, panel = panel.fun, 
					#ylim = c(-ylim.worm, ylim.worm), xlim = c(-xlim.worm, 
					#		xlim.worm), ylab = "Deviation", xlab = "Unit normal quantile", 
					show.given = show.given, bg = "wheat", pch = pch, 
					cex = cex, bar.bg = c(num = "light blue"), ...)
		#	if (overlap==0)
		#	{
		#		if  (total.points!=length(resid(object)))
		#			warning("the total number of points in the plot is not equal \n to the number of observations in y \n")
		#	}
		}
}
