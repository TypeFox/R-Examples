writecfg <- function(cfgfile = "autoquasimc.cfg",
					 verbose=0,
					 nruns=1,
					 samplingmethod="korobov",
					 log2nrays=20,
					 maxdepth=5,
					 russianroulette=c(0,0),
					 gridsize=c(10,10,10),
					 removeobjects="yes",
					 raysfromobjects="no",
					 onerayperspec="yes",
					 ignoredirect="no",
					 returnvar="yes",
					 outputfile="qmcarea.out",
					 returntype=c("F","D","U","L","H"),
					 lightmodel="Lambertian",
					 spectrumsamples=1,
					 sourcespectrum=c(660,1.0,730,1.0),
					 leaftop=c(0.1,0,0.1,0,1),
					 leafbottom=c(0.1,0,0.1,0,1),
				  	 sunazimuth = 0,
	                 sunaltitude = 90,
					 lightsourcefile="" 

){

	nrays <- 2^log2nrays
	returntype <- match.arg(returntype)
	
	if(is.na(outputfile))
		outputfile <- "no"
	else
		outputfile <- paste("yes",outputfile)
	
	
	# Solar direction (only one light source untik skyfile is more flexible).
	solarvec <- makesolarvec(sunaltitude,sunazimuth)

	
	vecp <- function(x,coll=" ")paste(x, collapse=coll)
	
	r <- c()
	r[1] <- paste("verbose:",verbose)
	r[2] <- paste("number of runs:",nruns)
	r[3] <- paste("sampling method:",samplingmethod)
	r[4] <- paste("number of rays:", nrays)
	r[5] <- paste("maximum depth:", maxdepth)
	r[6] <- paste("Russian roulette:",vecp(russianroulette))
	r[7] <- paste("grid size:",vecp(gridsize,coll=" x "))
	r[8] <- paste("remove objects:",removeobjects)
	r[9] <- paste("rays from objects:",raysfromobjects)
	r[10] <- paste("one ray per spectrum:",onerayperspec)
	r[11] <- paste("ignore direct light:",ignoredirect)
	r[12] <- paste("return variance:",returnvar)
	r[13] <- paste("output file:",outputfile)
	r[14] <- paste("return type:",returntype)
	r[15] <- paste("local light model:",lightmodel)
	r[16] <- paste("spectrum samples:",spectrumsamples)
	r[17] <- paste("source spectrum:",vecp(sourcespectrum))
	r[18] <- paste("leaf material (top):",vecp(leaftop))
	r[19] <- paste("leaf material (bottom):",vecp(leafbottom))
	
	if(lightsourcefile=="")
		r[20] <- paste("light source:",vecp(solarvec))
	else
		r[20] <- paste("light source file: ",lightsourcefile, sep="")
		
	writeLines(r, cfgfile)
	
return(invisible(r))
}















