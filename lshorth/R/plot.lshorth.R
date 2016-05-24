#$Id: plot.lshorth.R 92 2012-03-23 23:47:12Z gsawitzki $
plot.lshorth <- function(x, y, xlim = NULL, ylim = NULL,
		probs=NULL, 
		main="Shorth", 
		xlab=NULL, 
		ylab=NULL,
		frame.plot=TRUE,
		legendpos="topright", 
		rug=TRUE,
		rescale="neg", ...)
{
	stopifnot(inherits(x,"lshorth"))
	if (missing(ylim))  ylim <-NULL
	lshorthx <- x #fix me
	probs <- lshorthx$probs
	shorthm <- t(lshorthx$lshorth)
	if (is.null(xlab)) {
		if (!is.null(x$xname)){
            xlab <- paste(x$xname,", n=",length(lshorthx$x), sep="")
        } else{xlab <- paste(deparse(substitute(x), 50), 
        	", n=",length(lshorthx$x), collapse = "\n", sep="")}
        
    }
	if (is.null(rescale)) rsc <- "none" else{
	rsc<-match.arg(tolower(rescale),c("none","std","inv","neg"))
	}
	if (rsc=="std"){
		
	    shorthl <-min (lshorth(x=lshorthx$x,0.5,plot=FALSE)$lshorth)
	    # the length of the shorth
	    shorthmy <- shorthm/shorthl
	    
	    if (is.null(ylab)) {
		    ylab <- "std lshorth"
		    if(is.null(ylim)){
	ylim <- c(max(shorthmy)*1.1,0)
	if (is.na(ylim[2])) {
		ylim[2]<-0
	}}
	    }
	 }
     if (rsc=="inv") {
        shorthmy<- 1/shorthm
		if (is.null(ylab)) {
		    ylab <- "1/lshorth"
	    }
	    if (is.null(ylim)) {ylim<-1.1*range(shorthmy,finite=TRUE);ylim[1]<-0}
	 }
	 if (rsc=="neg") {
        shorthmy<- shorthm
		if (is.null(ylab)) {
		    ylab <- "lshorth"
	    }
	    if (is.null(ylim)){ylim<-c(1.1*range(shorthmy,finite=TRUE)[2],0)}
	 }

     if (rsc=="none") {
     	shorthmy<- shorthm
		if (is.null(ylab)) {
		    ylab <- "lshorth"
	    }
	    if (is.null(ylim)) {ylim<-range(shorthmy,finite=TRUE);ylim[1]<-0}
    }

	if (is.null(ylab)) {
		ylab <- "lshorth"
	}

	if (is.null(ylim)) ylim<-range(shorthmy,finite=TRUE)
	if (is.null(xlim)) xlim<-range(lshorthx$x[is.finite(lshorthx$x)])
	
	plot.new()
	plot.window(xlim=xlim,ylim=ylim, ...)
	axis(1)
	axis(2)
	title(main=main, xlab=xlab, ylab=ylab)
	if (frame.plot)box(...)
	if (rug) rug(lshorthx$x)
	
	lwd <- ceiling(6*(0.5-abs(probs-0.5)))
	
	for (px in 1:length(probs)){
	#	lwd <- ceiling(6*(0.5-abs(probs[px]-0.5)))
		lines(lshorthx$x,shorthmy[px,],lwd=lwd[px],...)
	}
	if (!is.null(legendpos)){
	temp <- legend(legendpos, legend = rep(" ",length(probs)),
               text.width = strwidth("0.0000"),
               lwd = lwd, xjust = 1, yjust = 1,
               title = expression(Coverage *" "* alpha),
               inset=0.05)
	text(temp$rect$left + temp$rect$w, temp$text$y,
     format(probs,digits=3), pos=2)
	}
	invisible(shorthm)
}

legend.lshorth <- function(legendpos,probs, ...){
	lwd <- ceiling(6*(0.5-abs(probs-0.5)))
	
	temp <- legend(legendpos, legend = rep(" ",length(probs)),
               text.width = strwidth("0.0000"),
               lwd = lwd, xjust = 1, yjust = 1,
               title = expression(Coverage *" "* alpha),
               inset=0.05, ...)
	text(temp$rect$left + temp$rect$w, temp$text$y,
     format(probs,digits=3), pos=2)
	}