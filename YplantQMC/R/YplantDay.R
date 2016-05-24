#'Run a simulation over a day with YplantQMC
#'
#'@description Interface to daily simulations with YplantQMC. Two objects are required to
#'run the simulation: a \code{plant3d} object, containing the plant structure
#'information, and a \code{met} object, containing weather data, solar
#'position, and number of timesteps.
#'
#'Optionally, a \code{phy} object is used which contains the leaf gas exchange
#'model for the simulation, to calculate photosynthesis (and possibly
#'transpiration rate) from light capture and other weather variables.
#'
#'Also optional is the use of a \code{hemi} object, which specifies shading by
#'a canopy.
#'
#'If you don't know where to start, run the example at the bottom of this page.
#'
#'See the arguments list above for the functions that are used to generate each
#'of the four objects. Note that the \code{plant} and \code{met} objects are
#'required, and \code{phy} and \code{hemi} are optional.
#'
#'This function is a user-friendly wrapper for \code{\link{runYplant}}. That
#'function should be used for all advanced simulations.
#'
#'@aliases YplantDay YplantDay.plant3d YplantDay.plant3dlist
#'@param x An object of class 'plant3d' or 'plant3dlist' (see
#'\code{\link{constructplant}} and \code{\link{readplantlist}}).
#'@param met An object of class 'ypmet', see \code{\link{setMet}}
#'@param phy An object of class 'ypphy', see \code{\link{setPhy}}
#'@param hemi An object of class 'yphemi', see \code{\link{setHemi}}
#'@param quiet If TRUE, does not write messages to the console.
#'@param writePSR If TRUE, writes a PSR output file.
#'@param PSRsuffix A suffix to be added to the PSR files that are written to disk.
#'@param writeOUT If TRUE, writes an OUT output file.
#'@param \dots Further arguments passed to \code{\link{runYplant}}
#'@return The \code{YplantDay} functions returns a list of class
#'\code{yplantsim}, which has \code{print} and \code{plot} methods (see
#'Examples).
#'
#'The list has the following components:
#'
#'\describe{ 
#'\item{plant}{The plant object used in the simulation}
#'\item{phy}{If provided, the phy object used in the simulation}
#'\item{hemi}{If provided, the hemi object used in the simulation}
#'\item{outdata}{A very lengthy dataframe with all results (see below)}
#'\item{nsteps}{Number of timesteps} 
#'\item{psrdata}{Totals and
#'averages by timestep (dataframe), see \code{\link{psrdata}}}
#'\item{met}{The met object used in the simulation} }
#'
#'The \code{outdata} dataframe in the \code{yplantsim} object lists results for
#'individual leaves, has the following variables.
#'
#'\describe{ 
#'\item{timeofday}{Time of day for current timestep (hours)}
#'\item{leafnr}{Leaf number} 
#'\item{timestep}{Length of current timestep (seconds)} 
#'\item{PAR0}{Above-canopy PAR}
#'\item{PARleaf}{Total PAR absorption} 
#'\item{PARdir}{Direct solar radiation PAR absorption} 
#'\item{PARdiff}{Diffuse PAR absorption} 
#'\item{reldiff}{Relative diffuse radiation absorption (0-1).} 
#'\item{reldir}{Relative direct radiation absorption (0-1).}
#'\item{LA}{Individual leaf area (mm2)} 
#'\item{LAproj}{Projected leaf area (mm2)} 
#'\item{LAsunlit}{Sunlit, or 'displayed' leaf area (mm2)} 
#'\item{A}{CO2 assimilation rate (mu mol m-2 s-1)}
#'\item{E}{Transpiration rate (mmol m-2 s-1)}
#'\item{gs}{Stomatal conductance (mol m-2 s-1)} 
#'\item{A0}{CO2 assimilation rate for a horizontal unshaded leaf (mu mol m-2 s-1)} 
#'}
#'
#'Where PAR is photosynthetically active radiation (mu mol m-2 s-1).
#'
#'The absorptions \code{reldiff} and \code{reldir} are relative to an unshaded
#'horizontal surface.
#'
#'To extract relative diffuse radiation absorption from an \code{yplantsim}
#'object, for example: 
#'\preformatted{ 
#'mysim <- YplantDay(myplant, mymet)
#'reldif<- mysim$outdata$reldiff 
#'}
#'@author Remko Duursma
#'@seealso \code{\link{runYplant}},\code{\link{ypreport}}
#'@references See \url{http://www.remkoduursma/yplantqmc}
#'@keywords misc
#'@examples
#'
#'
#'\dontrun{
#'# Set location,
#'southernfrance <- setLocation(lat=44)
#'
#'# A daily weather object, use a constant beam fraction of 0.4.
#'sunnyday <- setMet(southernfrance, month=6, day=21, nsteps=12, Tmin=9, Tmax=29, PARday=22,
#'	fbeamday=0.4, fbeammethod="constant")
#'
#'# Light response curve:
#'toonalrc <- setPhy("lightresponse", 
#'	leafpars=list(Amax=14.5, Rd=1.4, phi=0.05, theta=0.5, reflec=0.1, transmit=0.05))
#'
#'# Run YplantQMC for a day. Use the built-in 'largegap' hemiphoto.
#'toonarun <- YplantDay(toona, sunnyday, toonalrc, largegap)
#'}
#'
#'
#'@rdname YplantDay
#'
#'@export YplantDay
YplantDay <- function(x,...)UseMethod("YplantDay")

#'@rdname YplantDay
#'
#'@method YplantDay plant3dlist
#'@S3method YplantDay plant3dlist
YplantDay.plant3dlist <- function(x, met, phy=NULL, hemi=NULL, 
                                  ...){
  
  plants <- x
  nplants <- attributes(x)$nplants
  
  ypres <- list()
  for(i in 1:nplants){
    
    cat(paste(c(rep("=",35),"\n"),collapse=""))
    cat("Plant",i,"out of",nplants,"\n")
    cat(paste(c(rep("-",35),"\n"),collapse=""))
    ypres[[i]] <- YplantDay(plants[[i]], met, phy, hemi,  ...)
    
  }
  
  class(ypres) <- "yplantsimlist"
  
  return(invisible(ypres))
}


#'@rdname YplantDay
#'
#'@method YplantDay stand3d
#'@S3method YplantDay stand3d
YplantDay.stand3d <- function(x,...){
  
  stop("No method exists yet for a stand3d object. Will be added soon.")
  
}



#'@rdname YplantDay
#'
#'@method YplantDay plant3d
#'@S3method YplantDay plant3d
YplantDay.plant3d <- function(x, met, phy=NULL, hemi=NULL, quiet=FALSE,
                              writePSR=TRUE, PSRsuffix="", writeOUT=FALSE, ...){  
  
  plant <- x
  
  # Windows only. # MC 4/12/2012 - updated to include Mac OS X
  if(.Platform$OS.type != "windows" && (Sys.info()[['sysname']] != "Darwin"))
    stop("QuasiMC runs on Windows and Mac OS X only.")
  
  # Check input
  if(!inherits(met, "ypmet"))stop("Need object of class 'ypmet'; see ?setMet.")
  if(!is.null(phy) && !inherits(phy, "ypphy"))stop("Need object of class 'ypphy'; see ?setPhy.")
  if(!is.null(hemi) && !inherits(hemi, "yphemi"))stop("Need object of class 'yphemi'; see ?setHemi.")
  
  # Extract some data from met object.
  nsteps <- nrow(met$dat)
  daylength <- met$daylength
  hours <- met$dat$timeofday
  
  #
  if(!quiet){
    message("Yplant - simulation started.\n")
    cat(paste(c(rep("-",30),"\n"),collapse=""))  # a horizontal line.
    flush.console()
  }
  starttime <- proc.time()[3]
  
  # Run ray-tracer for every timestep.
  res <- list()
  message("Running QuasiMC and leaf gas exchange model.");flush.console()
  
  # First timestep: Calculate everything.
  res[[1]] <- runYplant(plant, phy, hemi, 
                        altitude=met$dat$altitude[1],
                        azimuth=met$dat$azimuth[1],PAR0=met$dat$PAR[1], 
                        Tair=met$dat$Tair[1], 
                        VPD=met$dat$VPD[1],Ca=met$dat$Ca[1], fbeam=met$dat$fbeam[1],
                        rewriteplantfile=TRUE,   # do not rewrite plant file; unless 1st timestep.
                        delfiles=FALSE,      # do not delete QuasiMC file; unless last timestep.
                        ...)
  othervars <- data.frame(timeofday = hours[1],
                          leafnr = plant$leafdata$leafnr,
                          timestep = ifelse(is.numeric(daylength),60*60 * daylength / nsteps,NA))
  
  res[[1]] <- cbind(othervars, res[[1]])
  reldiffPAR <- res[[1]]$reldiff
  if(!quiet){
    message("Timestep ",1," out of ",nsteps," completed.")
    flush.console()
  }
  
  # All other timesteps : only run QuasiMC if solar angle has changed (to save computing speed).
  if(nsteps > 1){
    for(i in 2:nsteps){
      
      # If solar altitude&azimuth same as in last timestep, use reldir from last timestep.
      if(met$dat$altitude[i] == met$dat$altitude[i-1] && 
           met$dat$azimuth[i] == met$dat$azimuth[i-1]){
        
        reldir <- res[[i-1]]$reldir
      } else {
        reldir <- NULL
      }
      
      res[[i]] <- runYplant(plant, phy, hemi, 
                            reldiff=reldiffPAR, 
                            reldir=reldir,
                            altitude=met$dat$altitude[i],
                            azimuth=met$dat$azimuth[i],PAR0=met$dat$PAR[i], 
                            Tair=met$dat$Tair[i], 
                            VPD=met$dat$VPD[i],Ca=met$dat$Ca[i], fbeam=met$dat$fbeam[i],
                            rewriteplantfile=FALSE,
                            delfiles=ifelse(i==nsteps, TRUE, FALSE),      # do not delete QuasiMC file; unless last timestep.
                            ...)
      if(!quiet){
        message("Timestep ",i," out of ",nsteps," completed.")
        flush.console()
      }
      
      othervars <- data.frame(timeofday = hours[i],
                              leafnr = plant$leafdata$leafnr,
                              timestep = ifelse(is.numeric(daylength),60*60 * daylength / nsteps,NA))
      
      res[[i]] <- cbind(othervars, res[[i]])
    }
  }
  
  # Return; should also include simulation settings in here.
  l <- list()
  l$plant <- plant
  if(!is.null(phy))l$phy <- phy
  if(!is.null(hemi))l$hemi <- hemi
  l$outdata <- do.call("rbind", res)
  l$nsteps <- nsteps
  
  # 'Plant Summary Report' data : sums and averages by timestep.
  makepsrdata <- function(out, plant, 
                          sumvars=c("A","A0","E","PARleaf","PAR0","PARinc","PARdiff","PARdir")){
    
    LAplant <- plant$leafarea
    
    # All are in (mu/m)mol m-2(leaf) s-1.
    sumvars <- intersect(names(out), sumvars)
    dfr1 <- out[,sumvars] * 10^-6 * out$LA
    dfr1$timeofday <- out$timeofday
    psrdata <- aggregate(dfr1,by=list(dfr1$timeofday),FUN=sum)  # mu mol s-1
    psrdata <- psrdata / LAplant  # convert back to mu mol m-2 s-1
    psrdata$timeofday <- as.vector(with(out, tapply(timeofday, timeofday, unique)))
    psrdata$timestep <- as.vector(with(out, tapply(timestep, timeofday, mean)))
    psrdata$LAplant <- LAplant
    
    # Projected and sunlit (displayed) leaf area.
    agla <- aggregate(out[,c("LAproj","LAsunlit")],by=list(out$timeofday),FUN=sum)*10^-6
    psrdata$LAproj <- agla[,2]
    psrdata$LAsunlit <- agla[,3]
    
    # clean up
    # Lame way to do this .... but it works, for now.
    tod <- psrdata$timeofday
    psrdata$timeofday <- NULL
    psrdata$Group.1 <- tod
    names(psrdata)[1] <- "timeofday"
    psrdata$Group.1 <- NULL
    
    return(psrdata)
  }
  l$psrdata <- makepsrdata(l$outdata, plant)
  
  
  l$met <- met
  endtime <- proc.time()[3]
  l$elapsedtime <- endtime - starttime 
  if(!quiet)
    message("\nSimulation completed in ", round(endtime-starttime,2)," seconds.")
  
  class(l) <- "yplantsim"
  
  if(writePSR)writePSRfile(l, suffix=PSRsuffix)
  if(writeOUT)writeOUTfile(l)
  
  
  return(l)
  
}

