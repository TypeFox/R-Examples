# Purposely hidden function; this is a developer utility.
checkHemi <- function(hemi){
	
  tur <- YplantQMC::turtle482
  gapfracs <- evalHemi(hemi, altitude=tur$altitude, azimuth=tur$azimuth,degrees=FALSE)	
	hX <- cos(tur$altitude) * sin(tur$azimuth)
	hY <- cos(tur$altitude) * cos(tur$azimuth)
	plot(hemi)
	points(hX,hY,col='lightgrey',cex=0.5, pch=19)
	points(hX,hY,col='red',cex=3*gapfracs$gapfraction)
}