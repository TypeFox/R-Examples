#'Summarize 3D plants
#'
#'Summarize the 3D plant in various useful ways. Requires an object of class
#'\code{plant3d}, made with \code{\link{constructplant}}.
#'
#'The \code{summary.plant3d} prints a number of plant summary variables. They
#'are also stored in a list. These are the variables that are currently
#'calculated:
#'
#'\preformatted{ 
#'LA - Total leaf area (m2) 
#'meanleafsize - Mean leaf size (cm2)
#'nleavesp - Number of leaves 
#'leaflen - Mean leaf length (cm) 
#'meanleafang - Mean leaf angle (deg) 
#'wmeanleafang - Mean leaf angle weighted by leaf area
#'(deg) Xellipsoid - Ellipsoidal leaf angle dist. par.  
#'crownvol - Crown volume (convex hull) (m3) 
#'crownsurf - Crown surface area (convex hull) (m2) 
#'ALAC - Crown density (AL/AC) (m2 m-2) 
#'cw - Crown width (m) 
#'cl - Crown length (m)
#'htot - Total height(m) 
#'cshape - Crown shape index (-) 
#'stemsurf - Stem + branch surface area (cm2) 
#'stemvol - Stem + branch volume (cm3) 
#'stemdiam - Stem base diameter (mm) 
#'meanpath - Mean pipe length (mm) 
#'sdpath - Standard deviation of pipe length (mm) 
#'totlen - Total woody segment length (mm) 
#'Ek - Expected distance to 5 nearest leaves (no edge corr.)  
#'Ek2 - Expected distance to 5 nearest leaves (with edge corr.)  
#'Ok - Observed distance to 5 nearest leaves 
#'disp - Dispersion parameter (no edge corr.)  
#'disp2 - Dispersion parameter (with edge corr.)  
#'STARbar - (Optional, only when calcSTARbar = TRUE).  }
#'
#'Note that when generating a summary on a \code{plant3dlist} object, the above
#'information is written to an outputfile. The outputfile is tab-delimited, and
#'has the name 'Plant summaries-YYYY-MM-DD.txt', using the current date.
#'
#'The following functions are called to calculate some of the summary
#'variables: \code{\link{leafdispersion}} for Ek,Ek2,Ok,disp,disp2;
#'\code{pathlen} (hidden function) for meanpath,sdpath,totlen; the
#'\code{\link{fitdistribution}} function in the \code{LeafAngle} package for
#'Xellipsoid; \code{\link{getR}} for cw; and \code{\link{crownhull}} for
#'crownsurf, crownvol and ALAC. The remainder of the variables are trivial
#'calculations from the plant input file (the .p file) or the leaf file.
#'
#'Optionally, the \eqn{\overline{STAR}} is calculated (when \code{calcSTARbar =
#'TRUE}), by calling \code{\link{STARbar}}.  See its help page for full
#'details. Note that all options for \code{\link{STARbar}} can also be set with
#'the summary function (see Examples).
#'
#'@aliases summary.plant3d summary.plant3dlist
#'@param object Object of class 'plant3d' (one 3D plant) or 'plant3dlist' (a
#'list of 3D plants).
#'@param nKErepeat Number of replicates for \code{\link{leafdispersion}}, see
#'its help page.
#'@param nsignif Number of digits for output (only for printing).
#'@param writefile If TRUE, writes a text file with the summary results in the
#'cur. working dir.
#'@param calcSTARbar If TRUE, also calculates STARbar and adds it to the
#'summary result.
#'@param \dots Further arguments passed to \code{\link{STARbar}}.
#'@return A list with components described above. See also the Examples below.
#'@author Remko Duursma
#'@seealso \code{\link{crownhull}}, \code{\link{leafdispersion}},
#'\code{\link{getR}}
#'@keywords misc
#'@examples
#'
#'
#'# Print summary (use built-in Toona plant):
#'summary(toona)
#'
#'# Or save summary as a list, access single values:
#'plantsumm <- summary(toona)
#'plantsumm$meanpath  # mean path length from soil to leaf
#'# See table above for names of single variables.
#'
#'\dontrun{
#'# Summary on a plant3dlist ('myplants' is constructed with 'readplantlist').
#'summary(myplants, writefile=TRUE)
#'
#'# Also calculate STARbar (with the exact method).
#'summary(myplant, calcSTARbar=TRUE, method="exact")
#'}
#'
#'@method summary plant3d
#'@S3method summary plant3d
#'@rdname summary.plant3d
#'@importFrom LeafAngle fitdistribution
summary.plant3d <- function(object, nKErepeat=10, nsignif=3, calcSTARbar=FALSE, ...){

	plant <- object
    nleaves <- plant$nleaves
	pdata <- plant$pdata
	if(plant$inputformat == "P")
		pdatL <- plant$pdata[plant$pdata$Lt >= 1, ]
	else # if(plant$inputformat == "Q")
		pdatL <- plant$qdata
	
	if(nleaves > 0){
	
	XYZ <- do.call("rbind", lapply(plant$leaves, function(x)x$XYZ))
	
	# Crown surface and volume
    if(nleaves > 3){
		ch <- crownhull(XYZ, plotit=FALSE)
		crownvol <- ch$crownvolume * 1E-09  # m3
		crownsurf <- ch$crownsurface * 1E-06  # m2
	} else {
		crownvol <- NA
		crownsurf <- NA
    }

	# Crown radius, length and shape factor
	meanR <- getR(XYZ)
	cw <- 2*meanR / 1000
	cl <- (max(XYZ[,3]) - min(XYZ[,3])) / 1000
	cshape <- crownsurf / (cw*cl)
	
	# total plant height
	htot <- max(XYZ[,3])/1000
	
	# leaf dispersion
	ld <- leafdispersion(plant, crownvol = crownvol, nleaves=nleaves)
	if(all(is.na(ld))){
		disp <- disp2 <- Ek <- Ek2 <- Ok <- NA
	} else {
		disp <- ld$disp_noedge
		disp2 <- ld$disp_edge
		Ek <- ld$Ek_noedge
		Ek2 <- ld$Ek_edge
		Ok <- ld$Ok
	}
	
	# mean leaf angle
  angs <- getangles(plant)
	meanang <- mean(angs)
	
	# leaf angle weighed by leaf area (more appropriate to relate to PA, for example)
	L <- pdatL$L.3
	wmeanang <- weighted.mean(angs, L^2)
    
	# Fit ellipsoidal distribution
	if(nleaves > 3){
		X <- fitdistribution(angs, "ellipsoid")$distpars
	} else {
		X <- NA
  }
    
	# Total leaf area:
	LA <- plant$leafarea
	meanleafsize <- 10^4 * LA / nleaves  # cm
	leaflen <- 0.1 * mean(L)
	
	# stem surface area (m2), stem volume (m3),
	# stem base diameter (mm), mean (and SD) path length from leaf to stembase (mm).  
	# and path length from leaf to soil.
	if(plant$inputformat == "P"){
		stemsurf <- 10^-6 * with(pdata, sum(D * pi * L) + sum(D.2 * pi * L.2))
		stemvol <- 10^-9 * with(pdata, sum((D/2)^2 * pi *L) + sum((D.2/2)^2 * pi * L.2))
		stembasediam <- pdata$D[1]
		if(stembasediam == 0)stembasediam <- pdata$D[2]
		if(nleaves >0){
			path <- try(pathlen(plant), silent=TRUE)
			if(inherits(path, "try-error")){
				path <- data.frame(totlen=NA)
				warning("Could not calculate path length for ",plant$pfile)
			}
		} else { 
			path <- data.frame(totlen=NA)
		}
		meanpath <- mean(path$totlen,na.rm=TRUE)
		sdpath <- sd(path$totlen,na.rm=TRUE)
		totlen <- sum(pdata$L) + sum(pdata$L.1)
	} else {
		stemsurf <- NA
		stemvol <- NA
		stembasediam <- NA
		meanpath <- NA
		totlen <- NA
		sdpath <- NA
	}	
	
	
	# For plants without leaves (yes, that's right).
	} else {
	
		stemsurf <- 10^-6 * with(pdata, sum(D * pi * L) + sum(D.2 * pi * L.2))
		stemvol <- 10^-9 * with(pdata, sum((D/2)^2 * pi *L) + sum((D.2/2)^2 * pi * L.2))
		stembasediam <- pdata$D[1]
		
		crownvol <- crownsurf <- leaflen <- 
		meanang <- wmeanang <- X <- Ek <- Ek2 <- Ok <- disp <- disp2 <- 
		meanpath <- sdpath <- totlen <- cw <- cl <- cshape <- LA <- meanleafsize <-
		NA
	
	} # end if(nleaves > 0)
	
	# STARbar calculation, if requested.
	if(calcSTARbar){
		run <- STARbar(plant,...) 
		STARbar_est <- run$DAbar / run$LA
	} else STARbar_est <- NA
	
	# (convex) crown projected area.
	AP <- 10^-6 * Silhouette(plant,azimuth=0,altitude=90)$H
	
	LFILE <- if(is.character(plant$lfile))plant$lfile else "unknown"
	
	obj <- list(plant$pfile, LFILE, crownvol,crownsurf,AP,nleaves,leaflen,
		meanang,wmeanang,X,Ek,Ek2,Ok,disp,disp2,
		stemsurf,stemvol,stembasediam,meanpath,
		sdpath,totlen,cw,cl,htot,cshape,LA,meanleafsize,STARbar_est)
	names(obj) <- c("pfile","lfile","crownvol","crownsurf","crownproj","nleavesp","leaflen",
		"meanleafang","wmeanleafang",
		"Xellipsoid","Ek","Ek2","Ok","disp","disp2",
		"stemsurf","stemvol","stemdiam","meanpath","sdpath","totlen",
		"cw","cl","htot","cshape","LA","meanleafsize","STARbar")
		
	obj$nsignif <- nsignif
		
	class(obj) <- "summary.plant3d"
	
	return(obj)
}


