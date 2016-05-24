"plot.scan.gwaa.2D" <-
function(x,y,...,df=1) {
	if (df == "all") {
#  		mx <- floor(max(-log10(x$P1df),-log10(x$P2df),na.rm=TRUE))
		image(x=x$map,y=x$map,z=-log10(x$P1df),...)
		image(x=x$map,y=x$map,z=-log10(t(x$P2df)),add=TRUE,...)
	} else if (df==2) {
#  		mx <- floor(max(-log10(x$P2df),na.rm=TRUE))
#		image(x=x$map,y=x$map,z=-log10(x$P2df),col=heat.colors(mx))
		image(x=x$map,y=x$map,z=-log10(x$P2df),...)
	} else {
#  		mx <- floor(max(-log10(x$P1df),na.rm=TRUE))
#		image(x=x$map,y=x$map,z=-log10(x$P1df),col=heat.colors(mx))
		image(x=x$map,y=x$map,z=-log10(x$P1df),...)
	}
}
