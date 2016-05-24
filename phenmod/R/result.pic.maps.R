result.pic.maps <- function(values, picPath=getwd(), 
				picName="beech-budburst", 
				silent=FALSE, 
				createFile=TRUE){
	
	if (!silent) { cat("Creating Picture of Budburst-DoY's calculated by model..") }

	#bb-doy-pim
	if (createFile){ png(filename=paste(picPath,"/",picName,"-doy.png",sep=""), 
		width=2000, height=2000, res=100) }
	trellis.par.set(fontsize=list(text=24,points=12), 
		par.main.text=list(cex=1), 
		superpose.symbol=list(cex=rep(1,7)),
		plot.symbol=list(cex=1),
		dot.symbol=list(cex=1, pch=18),
		axis.text=list(cex=1))
	picture <- levelplot(values$doy.model ~ values$x + values$y, values, 
		at=c(-10000,0,95,100,105,110,115,120,125,130,135,366), contour=FALSE,
		col.regions=c("grey",brewer.pal(n=10, name="PuOr")),
		xlab="East (in m)",ylab="North (in m)", main="Budburst-DoY calculated by PIM",
		colorkey=list(at=c(0,20,30,40,50,60,70,80,90,100,120), col=brewer.pal(n=10, name="PuOr"),
			labels=list(labels=c(0,95,100,105,110,115,120,125,130,135,365),
			at=c(0,20,30,40,50,60,70,80,90,100,120)))
		)
	print(picture)
	if (createFile){ dev.off() }
	if (!silent) { cat(" Done!\n") }
	
	if (!silent) { cat("Creating Picture of observed Budburst-DoY's..") }
	#bb-doy-observed
	if (createFile){ png(filename=paste(picPath,"/",picName,"-doy-observed.png",sep=""), 
			width=2000, height=2000, res=100) }
	trellis.par.set(fontsize=list(text=24,points=12), 
		par.main.text=list(cex=1), 
		superpose.symbol=list(cex=rep(1,7)),
		plot.symbol=list(cex=1),
		dot.symbol=list(cex=1, pch=18),
		axis.text=list(cex=1))
	picture <- levelplot(values$doy.observed ~ values$x + values$y, values, 
		at=c(-10000,0,95,100,105,110,115,120,125,130,135,366), contour=FALSE,
		col.regions=c("grey",brewer.pal(n=10, name="PuOr")),
		xlab="East (in m)",ylab="North (in m)", main="Observed Budburst-DoY",
		colorkey=list(at=c(0,20,30,40,50,60,70,80,90,100,120), col=brewer.pal(n=10, name="PuOr"),
			labels=list(labels=c(0,95,100,105,110,115,120,125,130,135,365),
			at=c(0,20,30,40,50,60,70,80,90,100,120)))
		)
	print(picture)
	if (createFile){ dev.off() }
	if (!silent) { cat(" Done!\n") }

	
	if (!silent) { cat("Creating Picture of \"Difference to Observation\" ..") }
	#difference of bb-doy calculated by PIM to bb-doy observed
	if (createFile){ png(filename=paste(picPath,"/",picName,"-difToObservation.png",sep=""),
			width=2000, height=2000, res=100) }
	trellis.par.set(fontsize=list(text=24,points=12), 
		par.main.text=list(cex=1), 
		superpose.symbol=list(cex=rep(1,7)),
		plot.symbol=list(cex=1),
		dot.symbol=list(cex=1, pch=18),
		axis.text=list(cex=1))
	picture <- levelplot(values$doy.dif ~ values$x + values$y, values, 
		at=c(-10000,-100,-20,-15,-10,-5,0,5,10,15,20,100), contour=FALSE,
		col.regions=c("grey",brewer.pal(n=10, name="PuOr")),
		xlab="East (in m)",ylab="North (in m)", 
		main="Difference of Budburst-DoY calculated by model to Budburst-DoY observed",
		colorkey=list(at=c(-25,-20,-15,-10,-5,0,5,10,15,20,25), col=brewer.pal(n=10, name="PuOr"),
			labels=list(labels=c(-100,-20,-15,-10,-5,0,5,10,15,20,100),
			at=c(-25,-20,-15,-10,-5,0,5,10,15,20,25)))
		)
	print(picture)
	if (createFile){ dev.off() }
	if (!silent) { cat(" Done!\n") }

	return(TRUE)
}
