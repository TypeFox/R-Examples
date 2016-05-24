#Author: Antonio, Fabio Di Narzo. Last Modified $Date: 2005-12-02 17:15:39 +0100 (ven, 02 dic 2005) $
stplot <- function(series, m, d, idt=1, mdt) {
	checkEmbParms(series, m, d)
	eps.max <- diff(range(series))*sqrt(m)
	res <- matrix(0, 10, mdt)
	res <- .C("stplot", series=as.double(series), length=as.integer(length(series)), m=as.integer(m), d=as.integer(d), mdt=as.integer(mdt), idt=as.integer(idt),eps.max=as.double(eps.max), res=as.double(res), PACKAGE="tseriesChaos")[["res"]]
	stp <- matrix(res, 10, mdt)
	eps.m <- min(stp)
	eps.M <- max(stp)
	plot(0, xlim=c(0, mdt*idt/frequency(series)), ylim=c(eps.m*0.99, eps.M*1.01), xlab="time", ylab="distance", type="n", main="Space-time separation plot")
	x <- seq(1/frequency(series), mdt*idt/frequency(series), by=idt/frequency(series))
	for(i in 1:10) lines(x,stp[i,])
	invisible(stp)
}
