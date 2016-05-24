runQMCUOC <- function(plant, hemi=NULL, reflec=NULL, transmit=NULL, 
	...){

		# No reflection or transmission set: black leaves.
		if(is.null(reflec) & is.null(plant$phy))reflec <- c(0,0)
		if(is.null(transmit) & is.null(plant$phy))transmit <- c(0,0)
		
		# Read from phy object, if available.
		if(is.null(reflec) & !is.null(plant$phy) & "reflec" %in% names(plant$phy$leafpars))
			reflec <- plant$phy$leafpars$reflec
		if(is.null(transmit) & !is.null(plant$phy) & "transmit" %in% names(plant$phy$leafpars))
			transmit <- plant$phy$leafpars$transmit

		if(length(reflec) == 1)reflec <- rep(reflec,2)
		if(length(transmit) == 1)transmit <- rep(transmit,2)
			
		# Turtle diffuse sky (uniformly distributed points on the sky hemisphere).
		writeturtlediffuse(hemi)
		diffile <- "turtlediffuse.dat"
		# diffile <- "sky_uniform.dat"
		
		# Config file for QuasiMC
		writecfg(cfgfile="autoquasimc.cfg", 
				outputfile="qmcarea.out",
				returntype="D",
				leaftop=c(reflec[1],0,transmit[1],0,1),
				leafbottom=c(reflec[2],0,transmit[2],0,1),
				lightsourcefile=diffile, ...)
			
		# Extract QMC file from the plant object (it is generated in constructplant()).
		infile <- plant$qmcinputfile
		outfile <- plant$qmcoutputfile
		writeLines(plant$qmcinput, infile)
			
		# Run model (only once; averages over all viewing angles).
		runquasimc("autoquasimc.cfg",infile,outfile,debug=FALSE)

qmcarea <- read.table("qmcarea.out", header=TRUE)

return(qmcarea)

}