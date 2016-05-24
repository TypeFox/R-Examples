plot.clustab <-
function(x,nmin=NULL,nmax=NULL, ...){
	if (!inherits(x, "clustab")) 
      	stop("use only with \"clustab\" objects")
	if (is.null(nmin)) nmin <- 2
	if (is.null(nmax)) nmax <- length(x$meanCR)+1
	meanCR <- x$meanCR[(nmin-1):(nmax-1)]
	plot(x=seq(nmin,nmax),meanCR,xaxt="n",ylim=c(0,1),xlab="number of clusters",ylab="mean adjusted Rand criterion",type="b",...)
	axis(side=1,at=seq(nmin,nmax),labels=paste(nmin:nmax))
}

