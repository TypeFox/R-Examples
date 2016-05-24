# File read_particles.R
# Part of the hydroPSO R package, http://www.rforge.net/hydroPSO/ ; 
#                                 http://cran.r-project.org/web/packages/hydroPSO
# Copyright 2011-2012 Mauricio Zambrano-Bigiarini & Rodrigo Rojas
# Distributed under GPL 2 or later

################################################################################
#                         'read_particles'                                     # 
################################################################################
# Author : Mauricio Zambrano-Bigiarini & Rodrigo Rojas                         #  
# Started: 08-Nov-2011,                                                        #
# Updates: 17-Feb-2012 ; 09-Mar-2012 ; 12-Nov-2012                             #        
################################################################################
# This function reads the 'Particles.txt' ouput file of hydroPSO and returns   #
# a data.frame with the position and fitness value of all the particles in the #
# swarm for all the iterations                                                 #
# IMPORTANT: the GoF of the particles HAS TO BE located in the 3rd column of   #
#            the file                                                          #
################################################################################

                       
read_particles <- function(file="Particles.txt", 
                           verbose=TRUE,
                           plot=TRUE,
                           # Plotting arguments
                           gof.name="GoF",
                           MinMax=NULL, 
                           beh.thr=NA, 
                           beh.col="red", 
                           beh.lty=1, 
                           beh.lwd=2, 
                           nrows="auto",
                           col="black",                            
                           ylab=gof.name, 
                           main=NULL,
                           pch=19, 
                           cex=0.5, 
                           cex.main=1.5,
                           cex.axis=1.5,
                           cex.lab=1.5,                      
                           #...,
                           breaks="Scott",
                           freq=TRUE,
                           do.pairs=FALSE, 
                           # Parameters for the 3D dotty plots
                           dp3D.names="auto", 
                           GOFcuts="auto", 
                           colorRamp= colorRampPalette(c("darkred", "red",
                            "orange", "yellow", "green", "darkgreen", "cyan")), 
                           alpha=1,
                           points.cex=0.7,
                           legend.pos="topleft", # not used yet,but included for avoiding warnings from 'plot_results'
                           #####################################################
                           #### PNG options
                           do.png=FALSE,
                           png.width=1500,
                           png.height=900,
                           png.res=90,
                           dotty.png.fname="Params_DottyPlots.png",
                           hist.png.fname="Params_Histograms.png",
                           bxp.png.fname="Params_Boxplots.png",
                           ecdf.png.fname="Params_ECDFs.png",
                           runs.png.fname="Params_ValuesPerRun.png",
                           dp3d.png.fname="Params_dp3d.png",
                           pairs.png.fname="Params_Pairs.png"
                           ) {                         
                         
  ##############################################################################
  # 1)                            Checkings                                    #
  ##############################################################################
  # Checking that 'file' exists
  if ( !file.exists(file) )
     stop( "Invalid argument value: The file '", basename(file), "' doesn't exist" )
     
  # Checking 'beh.thr'
  if ( !is.na(beh.thr) ) {
     if ( is.null(MinMax) )
        stop("Missing argument: 'MinMax' has to be provided before using 'beh.thr' !!")  
  } # IF end
  
  # Checking 'MinMax'
  if ( !is.null(MinMax) ) {
     if ( !(MinMax %in% c("min", "max")) )
        stop("Invalid argument: 'MinMax' must be in c('min', 'max')")
  } # IF end

  ##############################################################################
  # 2)                Reading ALL the PARAMETER SETS                           #
  ##############################################################################
  # Reading the file
  if (verbose) message( "                                                     ")  
  if (verbose) message( "[ Reading the file '", basename(file), "' ... ]" )  
  params <- read.table(file=file, header=TRUE, skip=0) 

  # computing the number of columns
  ncols <- ncol(params)  
  
  # Getting the goodness-of-fit of each particle
  gofs <- params[, 3]
  
  # Removing the first two colums with Iter and Part + the last column with GoF
  params <- params[, -c(1:3)] 

  # Amount of total parameter sets 
  nparamsets <- nrow(params)
  if (verbose) message( "[ Total number of parameter sets: ", nparamsets, " ]" )     

  # computing the number of parameters
  nparams <- ncol(params)
  
  # Filtering out those parameter sets above/below a certain threshold
  if (!is.na(beh.thr)) {  
     # Checking 'beh.thr'
     mx <- max(gofs, na.rm=TRUE)
      if (beh.thr > mx)
       stop("Invalid argument: 'beh.thr' must be lower than ", mx ,"!!")
    
    # Computing the row index of the behavioural parameter sets
    ifelse(MinMax=="min", beh.row.index <- which(gofs <= beh.thr), 
                          beh.row.index <- which(gofs >= beh.thr) )
    
    # Removing non-behavioural parameter sets and gofs
    params <- params[beh.row.index, ]
    gofs   <- gofs[beh.row.index]
   
    # Amount of behavioural parameter sets 
    nbeh <- nrow(params)
    if (verbose) message( "[ Number of behavioural parameter sets: ", nbeh, " ]" )
    
    # To avoid problems with 'plot_params'
    if (plot) beh.thr <- NA
  } # IF end 

  # Selecting best particle
  if ( !is.null(MinMax) ) {
    ifelse(MinMax=="min", best.row.index <- which.min(gofs), 
                          best.row.index <- which.max(gofs) )
    part.best <- params[best.row.index, ]
    gof.best  <- gofs[best.row.index]
  } else {
      part.best <- NA
      gof.best  <- NA
    } # ELSE end
  
  ##############################################################################
  # 3)                            Plotting                                     #
  ############################################################################## 
  
  if (plot) plot_particles(params=params, #parameter values
                           gofs=gofs,
                           gof.name=gof.name,
                           MinMax=MinMax, 
                           beh.thr=beh.thr, 
                           beh.col=beh.col, 
                           beh.lty=beh.lty, 
                           beh.lwd=beh.lwd, 
                           nrows=nrows,
                           col=col, 
                           ylab=ylab, 
                           main=main,
                           pch=pch, 
                           cex=cex, 
                           cex.main=cex.main,
                           cex.axis=cex.axis,
                           cex.lab=cex.lab,                     
                           #..., 
                           breaks=breaks,
                           freq=freq,
                           do.pairs=do.pairs,
                           # Parameters for the 3D dotty plots
                           dp3D.names=dp3D.names,
                           GOFcuts=GOFcuts,
                           colorRamp= colorRamp, 
                           alpha=alpha,
                           points.cex=points.cex,
                           legend.pos=legend.pos,
                           verbose=verbose,
                           #####################################################
                           #### PNG options
                           do.png=do.png,
                           png.width=png.width,
                           png.height=png.height,
                           png.res=png.res,
                           dotty.png.fname=dotty.png.fname,
                           hist.png.fname=hist.png.fname,
                           bxp.png.fname=bxp.png.fname,
                           ecdf.png.fname=ecdf.png.fname,
                           runs.png.fname=runs.png.fname,
                           dp3d.png.fname=dp3d.png.fname,
                           pairs.png.fname=pairs.png.fname
                           )  
  
  ##############################################################################
  # 4)                            Output                                       #
  ##############################################################################                   
  out <- list(part.params=params, part.gofs=gofs, best.param=part.best, best.gof=gof.best)
  
  return(out)
  
}  # 'read_particles' END



