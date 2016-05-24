# Write skyfile for QuasiMC; turtle UOC.
writeturtlediffuse <- function(hemi=NULL,npoints=482){
		
		if(npoints==482)
			sky <- YplantQMC::turtle482[,1:2]
		if(npoints==244)
# 			sky <- YplantQMC::turtle244[,1:2]
		if(npoints==59)
			sky <- YplantQMC::turtle[,1:2]
			
			
		if(is.null(hemi)){
			sky$intensity <- 1/nrow(sky)
		} else {
			sky$intensity <- evalHemi(hemi, altitude=sky$altitude, azimuth=sky$azimuth, 
				degrees=FALSE)$gapfraction
		}
		
		# QuasiMC expects zenith angle, not altitude.
		sky$altitude <- pi/2 - sky$altitude
		
		write.table(sky, "turtlediffuse.dat",col.names=FALSE, row.names=FALSE)
}
