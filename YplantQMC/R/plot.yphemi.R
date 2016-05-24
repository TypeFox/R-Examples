#'@method plot yphemi
#'@S3method plot yphemi
#'@rdname setHemi
plot.yphemi <- function(x,met=NULL,sungap=TRUE,
	projection=c("iso","flat"),warn=TRUE,bordercol='black', ...){

	projection <- match.arg(projection)
	hemi <- x
	
  o <- par(no.readonly=TRUE)
	on.exit(par(o))
  
  par(pty='s')
	plot(1, type='n',
	axes=FALSE, ann=FALSE,xlim=c(-1.15,1.15),ylim=c(-1.15,1.15))
	
	# Plot the tiles of the hemiphoto.
	for(i in 1:hemi$naz){
		for(j in 1:hemi$nalt){
			angs <- seq(hemi$azbins[i], hemi$azbins[i+1], length=25)
			if(projection == "flat"){
				r1 <- cos(hemi$altbins[j])
				r2 <- cos(hemi$altbins[j+1])
			}
			if(projection == "iso"){
				r1 <- 1 - hemi$altbins[j] / (pi/2)
				r2 <- 1 - hemi$altbins[j+1] / (pi/2)
			}
				
			x1 <- r1 * sin(angs)
			y1 <- r1 * cos(angs)
			x2 <- r2 * sin(angs)
			y2 <- r2 * cos(angs)
			xx <-c(x1,rev(x2))
			yy <-c(y1,rev(y2))
			
			polygon(xx,yy,col=grey(hemi$m[j, i]), border=bordercol)
		}
	}

	# Solar path.
	if(!is.null(met)){
		
		# Calculate solar az, alt, for many timesteps.
		hrs <- seq(met$sunrise+10/60, met$sunset-10/60, length=101)
		sunpos <- zenaz(met$year, met$month, met$day, 
			met$location$lat, met$location$long, met$location$tzlong,
			timeofday=hrs)
		sunpos2 <- zenaz(met$year, met$month, met$day, 
			met$location$lat, met$location$long, met$location$tzlong,
			timeofday=met$dat$timeofday)
		ALT <- sunpos$altitude
		AZ <- sunpos$azimuth
		
		if(projection == "flat"){
			sX <- cos(ALT*pi/180) * sin(AZ*pi/180)
			sY <- cos(ALT*pi/180) * cos(AZ*pi/180) 
			
			sunX <- cos(met$dat$altitude*pi/180) * sin(met$dat$azimuth*pi/180)
			sunY <- cos(met$dat$altitude*pi/180) * cos(met$dat$azimuth*pi/180) 
		}
		if(projection == "iso"){
			r <- 1 - ALT / 90
			sX <- r * sin(AZ*pi/180)
			sY <- r * cos(AZ*pi/180)
			
			rs <- 1 - met$dat$altitude / 90
			sunX <- rs * sin(met$dat$azimuth*pi/180)
			sunY <- rs * cos(met$dat$azimuth*pi/180)
		}
		
		gapfracdir <- evalHemi(hemi, met=met)$gapfraction
		
		if(sungap){
			points(sX, sY, col="darkorange2", pch=19, type='l')
			points(sunX, sunY, col="darkorange2", pch=19, cex=3*gapfracdir)
			if(max(gapfracdir,na.rm=TRUE) < 0.1 && warn)
				warning("Low gap fraction : solar path probably not visible. Try: sungap=FALSE.")
		} else {
			points(sX, sY, type='l', col="darkorange2")
		}

	}
	
	# Add circle; labels; legend.
	angs <- seq(0, 2*pi, length=101)
	x <- sin(angs)
	y <- cos(angs)
	points(x,y, type='l')
	text(0,1.1,expression(bold(N)))
	text(1.1,0,expression(bold(E)))
	text(0,-1.1,expression(bold(S)))
	text(-1.1,0,expression(bold(W)))
	fracs <- c(0.0,0.25,0.5,0.75, 1.0)
	legend("bottomleft", as.character(fracs), 
		fill=grey(fracs),cex=0.8,
		y.intersp=1.0, title="Gap fraction")
}
