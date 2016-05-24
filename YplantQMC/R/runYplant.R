#'A single simulation of YplantQMC
#'
#'@description Runs the YplantQMC model for one timestep. Runs the QuasiMC raytracer to
#'estimate absorbed PAR for every leaf on the plant, given diffuse and direct
#'radiation (set by \code{fbeam}, see below), the position of the sun, and
#'reflectance and transmittance of the foliage material.
#'
#'Required is a 3D plant object (see \code{\link{constructplant}}). Optionally,
#'a leaf gas exchange model is used (not needed if only light absorption is
#'calculated) to calculate photosynthesis (and optionally, transpiration rate).
#'Also optionally, a hemiphoto object is used to calculate shading by the
#'overstorey canopy.
#'
#'Output is not as easy to use as the more user-friendly
#'\code{\link{YplantDay}}. If you are only interested in diurnal simulations
#'(and plant totals by timestep), use that function. The \code{runYplant}
#'function is available for programming purposes (and more advanced custom
#'simulations).
#'
#'
#'The arguments \code{intern}, \code{debug}, \code{delfiles} and
#'\code{rewriteplantfile} should not be set by the user, unless you really know
#'what you are doing. These arguments exist for testing, and are used by
#'\code{YplantDay}.
#'
#'The arguments \code{reldiff} and \code{reldir} can be supplied if they are
#'already known (from a previous simulation, when the solar angle was the same,
#'in particular). If you are not sure, please do not set these arguments!
#'
#'@aliases runYplant runYplant.plant3d runYplant.stand3d
#'@param x An object of class 'plant3d', see \code{\link{constructplant}}
#'@param phy An object of class 'ypphy', see \code{\link{setPhy}}
#'@param hemi An object of class 'yphemi', see \code{\link{setHemi}}
#'@param reldiff Optional. A vector of relative diffuse absorption, same length
#'as number of leaves. See Details.
#'@param reldir Optional. A vector of relative direct absorption, same length
#'as number of leaves. See Details.
#'@param altitude,azimuth Solar altitude and azimuth (degrees).
#'@param fbeam Beam fraction (0-1). If 0, only diffuse interception is
#'calculated, if 1, only direct.
#'@param VPD Vapor pressure deficit (kPa)
#'@param PAR0 Incident PAR on a horizontal surface (mu mol m-2 s-1).
#'@param PARwhere If 'above', \code{PAR0} is given as an above-canopy value. If
#''below', it is below the canopy. See Details.
#'@param Ca Atmospheric CO2 concentration (ppm).
#'@param Tair Air temperature (deg C).
#'@param Patm Atmospheric pressure (kPa).
#'@param reflec Leaf reflectance (top, bottom of leaf).
#'@param transmit Leaf transmittance (top, bottom of leaf).
#'@param runphoto Whether to run leaf gas exchange model (default TRUE, or
#'FALSE when no phy object given).
#'@param intern If FALSE, returns output of QuasiMC to the console.
#'@param debug If TRUE, opens the QuasiMC debug window (for testing).
#'@param delfiles If TRUE, deletes intermediate files, and QuasiMC in/output
#'files.
#'@param rewriteplantfile If TRUE, writes the plant QuasiMC input file.
#'@param \dots Further arguments passed to \code{writecfg}.
#'@return This returns a dataframe with one row per leaf. The variables
#'included are (PAR is in units mu mol m-2 s-1):
#'
#'\describe{ \item{PAR0}{Incident PAR on a horizontal surface *above* the
#'canopy} \item{PARinc}{Incident PAR on a horizontal surface *below* the
#'canopy} \item{PARleaf}{Absorbed PAR (for each leaf)} \item{PARdir}{Absorbed
#'direct PAR} \item{PARdiff}{Absorbed diffuse PAR} \item{reldiff}{Relative
#'diffuse absorbed PAR (0 - 1)} \item{reldir}{Relative direct absorbed PAR (0 -
#'1)} \item{LA}{Leaf area (mm2)} \item{LAproj}{Projected leaf area (mm2)}
#'\item{LAsunlit}{Sunlit leaf area (mm2)} \item{A}{CO2 assimilation rate (mu
#'mol m-2 s-1)} \item{E}{Transpiration rate (mmol m-2 s-1)} \item{gs}{Stomatal
#'conductance (mol m-2 s-1)} \item{A0}{CO2 assimilation rate for a horizontal
#'leaf *below* the canopy.} }
#'@author Remko Duursma
#'@seealso \code{\link{YplantDay}}, \code{\link{setPhy}},
#'\code{\link{setHemi}}.
#'@references See \url{http://www.remkoduursma/yplantqmc}
#'@keywords misc
#'@examples
#'
#'
#'\dontrun{
#'
#'# Compare diffuse only to direct only
#'run_dir <- runYplant(pilularis, fbeam=1, altitude=90, azimuth=0, reflec=0.15, transmit=0.1)
#'run_diff <- runYplant(pilularis, fbeam=0, reflec=0.15, transmit=0.1)
#'
#'# Compare density functions of absorbed PAR by leaf:
#'plot(density(run_dir$PARleaf, from=0, to=1), xlim=c(0,1), main="", lwd=2, col="blue",
#'	xlab="Absorbed PAR (relative units)")
#'lines(density(run_diff$PARleaf, from=0, to=1), lwd=2, col="red")
#'legend("topright",c("Diffuse","Direct"), lwd=2, col=c("red","blue"))
#'}
#'
#' @export runYplant
#' @rdname runYplant
runYplant <- function(x,...)UseMethod("runYplant")

#'@method runYplant stand3d
#'@S3method runYplant stand3d
#'@rdname runYplant
runYplant.stand3d <- function(x,...){
  
  
  stand <- x  
  
  # List of leaves
  leaves <- list()
  for(i in 1:stand$nplants){
    leaves[[i]] <- stand$plants[[i]]$leaves
  }
  leaves <- do.call(c,leaves)
  
  # Fake 'plant' object; a link of all plants in the stand.
  p <- list(leaves=leaves)
  p$nleaves <- sum(stand$nleaves)
  
  # Make QuasiMC input file
  m <- makeQMCinput(p, writefile=FALSE)
  p$qmcinput <- m$qmcinput
  p$qmcinputfile <- m$inputfile
  p$qmcoutputfile <- m$outputfile
  
  p$phy <- NULL
  
  # $leafdata$area
  ld <- lapply(stand$plants, function(x)x$leafdata)
  ld <- do.call(rbind,ld)
  p$leafdata <- ld
  
  run <- runYplant.plant3d(p,...)
  
  run$plantnr <- rep(1:stand$nplants, stand$nleaves)
  
  return(run)
}


#'@method runYplant plant3d
#'@S3method runYplant plant3d
#'@rdname runYplant
runYplant.plant3d <- function(x, 
			phy=NULL,
			hemi=NULL,
			reldiff=NULL,  # vector of length nleaves (diffuse radiation, relative units).
			reldir=NULL,  # vector of length nleaves (direct radiation, relative units).
			# Solar parameters.
			altitude=90, 
			azimuth=0, 
			fbeam=1.0,
			# Environmental parameters.
			VPD=1.5, 
			PAR0=1,
			PARwhere=c("above","below"),
			Ca=390,
			Tair=25,
			Patm=101,
			# Leaf parameters.
			reflec=c(0.1,0.1),
			transmit=c(0.1,0.1),
			# Other options
			runphoto=TRUE,
			intern=TRUE,
			debug=FALSE,  
			delfiles=TRUE,
			rewriteplantfile=TRUE,   # Not for users! Only for YplantDay.
			...    # parameters to writecfg (for QuasiMC)
			){
	
	plant <- x
  
	# Windows only. # MC 7/11/2012 - updated to include Mac OS X
        # .Platform$OS.type returns "unix" for Mac, so need to check Sys.info()
	if((.Platform$OS.type != "windows") && (Sys.info()[['sysname']] != "Darwin"))
		stop("QuasiMC is currently available for Windows and Mac OS X only.")
	
	# If 'above', PAR0 is measured above the canopy (and so must be downscaled with hemi),
	# if 'below', PAR0 is already reduced by canopy shading.
	PARwhere <- match.arg(PARwhere)
	
	# Extract the phy object from the plant - if it exists.
	if(!is.null(plant$phy))phy <- plant$phy
	
	if(is.null(phy)){
		runphoto <- FALSE
	}
		
	# If the phy object has 'transmit' and 'reflec' : use those.
	if("reflec" %in% names(phy$leafpars))reflec <- phy$leafpars$reflec
	if("transmit" %in% names(phy$leafpars))transmit <- phy$leafpars$transmit	
	
	# Top and bottom reflec,transmit can be input.
	if(length(reflec) == 1)reflec <- rep(reflec,2)
	if(length(transmit) == 1)transmit <- rep(transmit,2)
	
	# Extract QMC file from the plant object (it is generated in constructplant()).
	infile <- plant$qmcinputfile		
	outfile <- plant$qmcoutputfile
	if(rewriteplantfile){
		writeLines(plant$qmcinput, infile)
	}
	
	# Input file for QuasiMC - DIRECT radiation only.
	writecfg(cfgfile="autoquasimc.cfg", 
		sunazimuth=azimuth, 
		sunaltitude=altitude, 
		outputfile="qmcarea.out",
		returntype="D",
		leaftop=c(reflec[1],0,transmit[1],0,1),
		leafbottom=c(reflec[2],0,transmit[2],0,1),
		...)
	# writecfg(cfgfile="autoquasimc.cfg", 
		# outputfile="qmcarea.out",
		# returntype="D",
		# leaftop=c(reflec[1],0,transmit[1],0,1),
		# leafbottom=c(reflec[2],0,transmit[2],0,1),
		# lightsourcefile="directqmcin.dat",
		# ...
	# )
		
	# Run the QuasiMC model for direct radiation (unless fbeam is zero - use only diffuse PAR).
	# Output is in relative units (horizontal unshaded area = 1.0).
	# Possibly, relative direct PAR absorption (reldir) was given as input,
	# based on output by this very function.
	if(is.null(reldir) && fbeam > 0){
		tm1 <- system.time(runquasimc("autoquasimc.cfg",infile,outfile,debug=debug,intern=intern))

		# Read output.
		# reldir <- readQMCout(plant$qmcoutputfile)  # E() module; obsolete.
		outres <- read.table("qmcarea.out", header=TRUE)
		reldir <- outres$rel_mean_flux0  # APAR relative to horizontal surface.
		LAproj <- outres$mean_projected_area
		LAsunlit <- outres$sunlit_area
		LAtot <- outres$actual_leaf_area
		
	} else {
		LAproj <- NA
		LAsunlit <- NA
		LAtot <- NA
		if(fbeam == 0)reldir <- 0
	}
	

	# Run QuasiMC for diffuse radiation (unless it is given as input).
	if(is.null(reldiff) && fbeam == 1.0){
		reldiff <- 0.0
	}
	if(is.null(reldiff) && fbeam < 1){
		# APAR relative to horizontal surface.
		reldiff <- runQMCUOC(plant, hemi=hemi, transmit=transmit, 
			reflec=reflec, ...)$rel_mean_flux0  
	}
	
	# Gap fraction.
	if(!is.null(hemi)){
		# Gap fraction for direct light
		gapfracdir <- evalHemi(hemi, altitude, azimuth)$gapfraction
		
		# Weighted mean gap fraction for diffuse light
		dif <- evalHemi(hemi, 
                    altitude=YplantQMC::turtle482$altitude, 
                    azimuth=YplantQMC::turtle482$azimuth,
			degrees=FALSE)
		gapfracdiff <- weighted.mean(dif$gapfraction, sin(dif$altitude))
		
		# Reduce PAR due to canopy shading (if PAR 'measured' above canopy).
		if(PARwhere == "above"){
			PARinc <- PAR0 * (fbeam*gapfracdir + (1-fbeam)*gapfracdiff)
		}
			
	} else {
		    
		PARinc <- PAR0
	
	}
		
	# Actual PPFD on leaves (vector of length nleaves).
	PARdir <- reldir * fbeam * PARinc
	PARdiff <- reldiff * (1 - fbeam) * PARinc
	PARs <- PARdir + PARdiff

	# Results by leaf.
	ypresults <- data.frame(PAR0=PAR0, PARinc=PARinc, PARleaf=PARs, PARdir=PARdir, 
		PARdiff=PARdiff, reldiff=reldiff, reldir=reldir)
	
  if(nrow(ypresults) != plant$nleaves)stop("Critical failure. nrow(yplant) != plant$nleaves")
  
	ypresults$LA <- plant$leafdata$area
	ypresults$LAproj <- LAproj
	ypresults$LAsunlit <- LAsunlit
	
	if(runphoto){
		# Run photosynthesis module.
		metvars <- list(PAR=PARs, Tair=Tair, Ca=Ca, VPD=VPD, Patm=Patm)
		photorun <- do.call(phy$leaffunction, c(phy$leafpars, metvars))

		# Run photosynthesis module for a flat unshaded leaf.
		metvars <- list(PAR=PARinc, Tair=Tair, Ca=Ca, VPD=VPD, Patm=Patm)
		photorunflat <- do.call(phy$leaffunction, c(phy$leafpars, metvars))

		# Add unshaded photosynthetic rate:
		photorun$A0 <- photorunflat$A
		
		# cbind.
		ypresults <- cbind(ypresults, photorun)
	} 
	
	if(delfiles){
		unlink(infile)
		unlink(outfile)
		unlink("qmcarea.out")
		unlink("autoquasimc.cfg")
		unlink("turtlediffuse.dat")
	}
	
return(ypresults)
}



	





	
