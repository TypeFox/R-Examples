dotplot.mtb <- function (x, xlim = NULL, main = NULL, xlab = NULL, ylab=NULL,
                     pch = 19, hist=FALSE, yaxis=FALSE, mtbstyle=TRUE)
{
    if (is.null(xlim))
        xlim <- range(pretty(range(x,na.rm=TRUE))) 
    if (is.null(main)) 
        main <- ""
    if (is.null(xlab)) 
        xlab <- ""
    if (is.null(ylab)) 
        ylab <- ""
    x <- sort(x)
    w <- table(x)
    if(hist) {
        x <- as.numeric(names(w))
        w <- unname(unclass(w))
    } else w <- unlist(lapply(w, function(n) { 1:n }))
    mw <- max(w)
    Nmax <- floor(par()$pin[2]/strheight("o",units="inches"))
    top <- if(mtbstyle) Nmax/2 else Nmax
    if(mw <= top & !hist) {
    	plot(range(x, na.rm = TRUE), c(0, 1), type = "n", xlab = "", 
             ylab = "", xlim = xlim, main = main, axes = FALSE)
	yr <- if(mtbstyle) c(-Nmax/2,Nmax/2) else c(0,Nmax)
    	par(usr = c(par()$usr[1:2], yr[1], yr[2]))
    	y <- strheight("o") * w
    	points(x, y, pch = pch)
    	axis(side = 1, pos = 0)
	if(xlab!="") axis(side=1,at=0.5*sum(xlim),pos=-2,labels=xlab,tick=FALSE)
	if(ylab!="") axis(side=2,at=0.5*yr[2],line=2,labels=ylab,tick=FALSE)
	if(yaxis) {
		b  <- max(1,ceiling(top/10))
		ll <- 1+floor(top/b)
		at <- seq(0,by=b,length=ll)
		axis(side = 2,at=at)
	}
    } else {
	nt <- mw+1
	yr <- if(mtbstyle) c(-nt,nt) else c(0,nt)
	if(hist) plot(x,w,type="h",xlab="",ylab="",xlim=xlim,
                      ylim=yr,main=main,axes=FALSE)
	else plot(x,w,pch=pch,xlab="",ylab="",xlim=xlim,
                  ylim=yr,main=main,axes=FALSE)
	pos <- if(mtbstyle) -0.02*mw else 0
	axis(side = 1,pos = pos)
	if(mtbstyle) {
		if(xlab!="") axis(side=1,at=0.5*sum(xlim),pos=-2,labels=xlab,tick=FALSE)
		if(ylab!="") axis(side=2,at=0.5*yr[2],line=2,labels=ylab,tick=FALSE)
	} else {
		if(xlab!="") title(xlab=xlab)
		if(ylab!="") title(ylab=ylab)
	}
	if(yaxis) {
		b  <- max(1,ceiling(nt/10))
		ll <- 1+floor(nt/b)
		at <- seq(0,by=b,length=ll)
		axis(side = 2,at=at)
	}
    }
}
