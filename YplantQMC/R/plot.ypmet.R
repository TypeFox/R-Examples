#'@method plot ypmet
#'@S3method plot ypmet
plot.ypmet <- function(x,...){

	dat <- x$dat
	
	plotcols <- c("darkorange2","red","forestgreen","darkgrey")
	
	par(yaxs="i", mar=c(5,7,4,7), mgp=c(1.6,0.6,0), tck=-0.01)
	with(dat, plot(timeofday, PAR, ylim=c(0,1.1*max(PAR)), 
		xlim=c(x$sunrise, x$sunset),
		xlab="Time of day (hours)",
		ylab=expression(PAR~~(mu*mol~m^-2~s^-1)),
		pch=19, type='o', col=plotcols[1], lwd=2))
	abline(v=x$sunrise, col="darkgrey", lty=5)
	abline(v=x$sunset, col="darkgrey", lty=5)
	box()
	par(new=TRUE)
	with(dat, plot(timeofday, Tair, col=plotcols[2],
		pch=19, type='o', axes=FALSE, ann=FALSE,lwd=2,
		ylim=c(0.9*min(Tair),max(Tair)*1.1),
		xlim=c(x$sunrise, x$sunset)))
	axis(2,line=3.5)
	mtext(expression(T[air]~~(""^"o"*C)),2,line=5.5)

	par(new=TRUE)
	with(dat, plot(timeofday, VPD, col=plotcols[3],
		pch=19, lwd=2, type='o', axes=FALSE, ann=FALSE,
		ylim=c(0,max(VPD)*1.1),
		xlim=c(x$sunrise, x$sunset)))
	axis(4)
	mtext("VPD (kPa)",4,line=1.6)
	
	par(new=TRUE)
	with(dat, plot(timeofday, fbeam, col=plotcols[4],
		pch=19, lwd=2, type='o', axes=FALSE, ann=FALSE,
		ylim=c(0.8*min(fbeam),max(fbeam)*1.4),
		xlim=c(x$sunrise, x$sunset)))
	axis(4, line=3.5)
	mtext("Beam fraction (-)",4,line=5.5)
	
	legend("bottom", c("PAR",expression(T[air]),"VPD","Beam frac."),
		col=plotcols, lty=1, lwd=2, pch=19)
	
	
}