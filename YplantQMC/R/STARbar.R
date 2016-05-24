#'Calculate displayed and projected leaf areas of a 3D plant
#'
#'@description From a 3D plant, calculates the area of all leaves facing the viewing angle
#'(from a given azimuth and altitude), or 'projected area', and the area of the
#'leaves that is visible ('displayed area'). By default, repeats calculations
#'from many viewing angles to estimate the hemi-spherically averaged STAR.
#'
#'This function calculates the displayed area (DA) and projected area (PA) of
#'the entire plant, either over some specified viewing angle(s), or (by
#'default) the hemispherically average values. From these averages, the
#'average STAR is estimated (see Duursma et al. 2012).
#'
#' The \code{\link{summary.plant3d}} function also calculates STARbar, if \code{calcSTARbar=TRUE} is set.
#'
#'There are four methods for the calculation of the displayed and projected
#'area from a given viewing angle: 
#'\describe{ 
#'\item{gridtracer}{A regular grid
#'is placed over the projected plant. Displayed area is calculated from the
#'number of intersecting grid points, and the grid size. The argument
#'\code{npside} sets the number of grid points on the shorter side of the
#'projected plant.  Slow for large plants, but good estimates are still
#'obtained with low values of \code{npside}.} 
#'\item{exact}{A polygon method to
#'calculate the exact displayed area. Useful for testing, but \emph{very}
#'slow.} 
#'\item{QuasiMC}{Runs \code{QuasiMC} over the entire hemisphere - the
#'\code{QuasiMC} model does the averaging.  Relatively slow for small plants,
#'but fast for large plants.} 
#'\item{slowquasimc}{Do not use this method (for
#'testing only). Runs \code{QuasiMC} once per viewing angle, then averages the
#'results.  Painfully slow.} }
#'
#'There are three integration methods. Note that if \code{method = "QuasiMC"},
#'the integration method is ignored. 
#'\describe{ 
#'\item{Turtlesky}{Integrates
#'over 58 points, that are placed approx. uniformly over the hemisphere.}
#'\item{Yplant}{Integrates over 160 angles, distributed over the hemisphere as
#'in the original Yplant.  Weighing is applied over these angles to yield the
#'average STAR.} 
#'\item{Turtle244}{As Turtlesky, but uses 244 points (slow method).} }
#'
#'@aliases STARbar STARbar.plant3dlist STARbar.plant3d raytrace
#'plot.tracedplant
#'@param object An object of class 'plant3d', see \code{\link{constructplant}}
#'@param method The method to calculate the displayed area. See Details.
#'@param integration The integration method to calculate the average STAR over
#'the hemisphere. See Details.
#'@param progressbar If TRUE (default), displays a graphical progress bar.
#'@param npside For \code{method = "gridtracer"}, the number of grid cells per
#'side for raytracing.
#'@param returnldr If TRUE, returns a dataframe with results per viewing angle.
#'@param quiet If TRUE, prints no progress to the screen.
#'@param silhouette If TRUE, also calculates the 2D convex hull around the
#'projected plant (see \code{\link{Silhouette}}).
#'@param azimuth,altitude Azimuth and altitude of the viewing direction
#'(Optional, in Radians).
#'@param \dots Further arguments are ignored, for now.
#'@author Remko Duursma
#'
#'
#'@seealso
#'\code{\link{projectplant}},\code{\link{constructplant}},\code{\link{Silhouette}},
#'\code{\link{summary.plant3d}}
#'@references Duursma, R.A., D.S. Falster, F. Valladares, F.J. Sterck, R.W.
#'Pearcy, C.H. Lusk, K.M. Sendall, M. Nordenstahl, N.C. Houter, B.J. Atwell, N.
#'Kelly, J.W.G. Kelly, M. Liberloo, D.T. Tissue, B.E. Medlyn and D.S.
#'Ellsworth. 2012.  Light interception efficiency explained by two simple
#'variables: a test using a diversity of small- to medium-sized woody plants.
#'New Phytologist. 193:397-408.
#'
#'See also \url{http://www.remkoduursma/yplantqmc}
#'@keywords misc
#'@examples
#'
#'\dontrun{
#'# Get STARbar for the built-in Toona plant:
#'# (Settings are fast, not accurate)
#'STARbar(toona, method="gridtracer", npside=15)
#'
#'# For exact STAR, use:
#'STARbar(toona, method="exact")
#'
#'# To produce an LDR file (as in the original Yplant), use:
#'clidstar <- STARbar(toona, method="gridtracer", npside=30, integration="Yplant", returnldr=T)
#'write.table(clidstar$ldr, "toona.LDR")
#'}
#'
#'
#'@rdname STARbar
#'@export STARbar
STARbar <- function(object,...){
	UseMethod("STARbar")
}

#'@rdname STARbar
#'
#'@method STARbar plant3dlist
#'@S3method STARbar plant3dlist
STARbar.plant3dlist <- function(object, quiet=FALSE, ...){
  
  plants <- object
  nplants <- attributes(plants)$nplants
  pfiles <- attributes(plants)$pfiles
  lfiles <- attributes(plants)$lfiles
  
  dapas <- list()
  for(i in 1:nplants){
    
    tm <- system.time(dapas[[i]] <- try(STARbar(plants[[i]],quiet=TRUE,...)))
    elaps <- unname(tm[3])
    
    if(inherits(dapas[[i]], "try-error")){
      warning("STARbar returned error for plant ",pfiles[i])
      dapas[[i]] <- NA
      next
    }
    if(!quiet && nplants>1){
      
      cat("Plant",i,"of",nplants,"done in", elaps ,"sec.\n")
      flush.console()
    }
    dapatmp <- as.data.frame(dapas[[i]][1:4]) 
    dapas[[i]]$elapsed <- elaps
    dapas[[i]]$pfile <- pfiles[i]
    dapas[[i]]$lfile <- lfiles[i]
  }
  
  if(length(dapas)==1)
    return(dapas[[1]])
  else {
    class(dapas) <- "STARbarlist"
    return(dapas)
  }
}


#'@rdname STARbar
#'
#'@method STARbar plant3d
#'@S3method STARbar plant3d
STARbar.plant3d <- function(object, 
                            method=c("gridtracer","exact","QuasiMC","slowquasimc"),
                            integration=c("Turtlesky","Yplant","Turtle244"), 
                            progressbar=TRUE, 
                            returnldr=FALSE, 
                            quiet=FALSE,
	                          npside=25, 
                            silhouette=TRUE, 
                            azimuth=NA, 
                            altitude=NA, 
                            ...
	){

	plant <- object
	method <- match.arg(method)
	integration <- match.arg(integration)

        
	if(method %in% c("QuasiMC","slowquasimc") && .Platform$OS.type != "windows" && 
       (Sys.info()[['sysname']] != "Darwin")){
	
    stop("Select different method: QuasiMC is currently only available in Windows and Mac OS X.")
	
	}
  
	if(quiet)progressbar <- FALSE
	
	# Get viewing angles to calculate DA and PA from.
	if(all(is.na(altitude)) | all(is.na(azimuth))){
		if(integration == "Yplant"){
			altaz <- YplantQMC::yplantaltaz
		}
		if(integration == "Turtlesky"){
			altaz <- YplantQMC::turtle
		}
		if(integration == "Turtle244"){
			altaz <- YplantQMC::turtle244
		}
		AZ <- altaz$azimuth
		ALT <- altaz$altitude
		nangles <- length(AZ)
		sphericalSTAR <- TRUE
	} else {
		AZ <- azimuth
		ALT <- altitude
		nangles <- length(AZ)
		quiet <- TRUE
		integration <- "Turtlesky"  # because STARbar is then simply average of specified angles.
		sphericalSTAR <- FALSE
		if(method=="QuasiMC")
			stop("Cannot use QuasiMC to calculate STAR for pre-defined angles (yet).")
		
		if(!quiet)message("Calculating STAR for pre-specified angles.")
	}

	# Total plant leaf area.
	LA <- plant$leafarea*10^6
	
	# init.
	DAs <- PAs <- Hs <- vector("numeric", nangles)
	
	# Msg.
	if(!quiet){
		message("Calculating STARbar with method: ",method,", integration scheme: ",integration,".")
		flush.console()
	}
	
  
	if(method == "QuasiMC"){ 
	
		# 160 angles not yet supported (will perhaps be?)
		if(integration=="Yplant")
      stop("QuasiMC with Yplant integration not implemented.")

		# Run QMC in UOC mode (with 'black' leaves).
		qmcdapa <- runQMCUOC(plant, transmit=0, reflec=0,...)
	
		# Total displayed & projected area, STARbar.
		DA <- sum(qmcdapa$sunlit_area)
	  PA <- sum(qmcdapa$mean_projected_area)
	  LA <- sum(qmcdapa$actual_leaf_area)
		starbar <- DA / LA
		pabar <- PA / LA
		l <- list(STARbar = starbar, DAbar = DA, PAbar = PA, LA = LA, PALAbar = pabar)
		l$ldr <- NA
		l$sphericalSTAR <- sphericalSTAR
	}
	
	
	# Run QuasiMC once for every angle. (Slow, for testing only).
	if(method == "slowquasimc"){
	
		# Initialize.
		if(progressbar){
			wp <- txtProgressBar(title = "Calculating STAR", 
			label = "", min = 0, max = nangles, initial = 0, width = 50,style=3)
		}	
	
		# Don't do Yplant integration - will be much too slow.
		if(integration %in% c("Yplant","Turtle244"))
			stop("slowquasimc with Yplant integration not implemented (too slow).")
		
		# Extract QMC file from the plant object (it is generated in constructplant()).
		infile <- plant$qmcinputfile
		outfile <- plant$qmcoutputfile
		writeLines(plant$qmcinput, infile)
	
		for(i in 1:nangles){
		
			res <- DAPAQMC(plant,AZ[i],ALT[i], infile,outfile)
			DAs[i] <- res$DA
			PAs[i] <- res$PA
			Hs[i] <- res$H
		
		if(progressbar)setTxtProgressBar(wp, i)
		}
		if(progressbar)close(wp)
		
		# The LDR format.
		ldr <- data.frame(azimuth=AZ, altitude=ALT, DA=DAs, PA=PAs)
		
		# If turtle sky, each grid point has the same weight. Simple averages.
		DAbar <- mean(DAs)
		PAbar <- mean(PAs)
		
		starbar <- DAbar / LA
		pabar <- PAbar / LA
		l <- list(STARbar = starbar, DAbar = DAbar, PAbar = PAbar, LA = LA, PALAbar = pabar)
		l$ldr <- ldr
		l$sphericalSTAR <- sphericalSTAR
	}
	
  if(method == "exact")
    stop("Method temporarily unavailable. Contact package maintainer.")
  
	if(method == "gridtracer"){
		
		# Initialize.
		if(progressbar){
			wp <- txtProgressBar(title = "Calculating STAR", 
			label = "", min = 0, max = nangles, initial = 0, width = 50,style=3)
		}	
 
		# Calculate DA, PA for all angles.
		for (i in 1:nangles) {
		
			res <- DAPA(plant, azimuth=AZ[i], altitude=ALT[i], exact=FALSE, npside=npside)
			DAs[i] <- res$DA
			PAs[i] <- res$PA
			Hs[i] <- res$H
			
		if(progressbar)setTxtProgressBar(wp, i)
		}
		if(progressbar)close(wp)
		
		# The result is identical to the LDR file in Yplant.
		ldr <- data.frame(azimuth=AZ, altitude=ALT, DA=DAs, PA=PAs, H=Hs)
		
		# Make average, Duursma et al. 2012 (NewPhyt) method.
		if(integration == "Yplant"){
			ldr$az <- rep(1:8,20)
			ldr$zen <- rep(1:20, each=8)
			zenmeans <- aggregate(ldr[,c("DA","PA")], list(ldr$zen), mean)
			s <- seq(0,90,length=21)
			zenmeans$alpha <- (pi/180)*(s[2:21] + s[1:20])/2
			
			da <- 2.25*pi/180
			zenmeans$Weight <- with(zenmeans, sin(alpha+da) - sin(alpha-da))
			
			DAbar <- with(zenmeans, weighted.mean(DA, Weight))
			PAbar <- with(zenmeans, weighted.mean(PA, Weight))
			ldr$az <- ldr$zen <- NULL
		}	
		
		# If turtle sky, each grid point has the same weight. Simple averages.
		if(integration %in%  c("Turtlesky","Turtle244")){
			DAbar <- mean(ldr$DA)
			PAbar <- mean(ldr$PA)
		}
		
	starbar <- DAbar / LA
	pabar <- PAbar / LA
	l <- list(STARbar = starbar, DAbar = DAbar, PAbar = PAbar, LA = LA, PALAbar = pabar)
	l$ldr <- ldr
	l$sphericalSTAR <- sphericalSTAR
	}

	class(l) <- "STAR"
return(l)
}

