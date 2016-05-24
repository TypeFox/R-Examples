# Displayed area (DA) and projected area (PA) using QuasiMC (QMC),
# for one viewing angle at a time.
DAPAQMC <- function(plant, azimuth, altitude, inputfile, outputfile){

	k <- 180/pi
	writecfg(cfgfile="autoquasimc.cfg", sunazimuth=k*azimuth, 
		sunaltitude=k*altitude, outputfile="qmcarea.out")
	
	runquasimc("autoquasimc.cfg",inputfile,outputfile,debug=FALSE)
	
	qmcdapa <- read.table("qmcarea.out", header=TRUE)

	l <- list(DA = sum(qmcdapa$sunlit_area),
	          PA = sum(qmcdapa$mean_projected_area),
			  azimuth=azimuth*k,
			  altitude=altitude*k,
			  exact=FALSE,
			  LA=sum(qmcdapa$actual_leaf_area))
	
	# 
	l$H <- Silhouette(plant, azimuth*k, altitude*k)$H
	
return(l)	

}
