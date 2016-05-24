`plot.apos` <-
function (x,label=TRUE,grid=TRUE,type="n",...) 
{
	if (inherits(x,"eqc")) { xlim=c(0,24); xaxp=c(0,24,6); gr=0:12*2;
		main="Equatorial Coordinates";
		xlab="Right Ascension";
		ylab="Declination";
		}
	else
	{ xlim=c(0,360); xaxp=c(0,360,6); gr=0:12*30; 
		if (inherits(x,"hoc")) {
			main="Horizontal Coordinates";
			xlab="Azimuth";
			ylab="Altitude";
			}
		if (inherits(x,"ecc")) {
			main="Ecliptic Coordinates";
			xlab="Longitude";
			ylab="Latitude";
			}
	}

	plot(as.vector(x[,1]),as.vector(x[,2]),type=type,
		main=main,
		xlab=xlab,
		ylab=ylab,
		xlim=xlim,
		xaxp=xaxp,
		yaxp=c(-90,90,12),
		ylim=c(-90,90),
		cex=0.75,...);
	if (grid) {abline(h=0); abline(v=gr,lty=3);}
	if (label) text(as.vector(x[,1]),as.vector(x[,2]),labels=rownames(x),cex=0.6);

}

