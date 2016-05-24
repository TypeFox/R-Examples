# File read_results.R
# Part of the hydroPSO R package, http://www.rforge.net/hydroPSO/ ; 
#                                 http://cran.r-project.org/web/packages/hydroPSO
# Copyright 2011-2012 Mauricio Zambrano-Bigiarini & Rodrigo Rojas
# Distributed under GPL 2 or later

################################################################################
#                         'read_results'                                       # 
################################################################################
# Author : Mauricio Zambrano-Bigiarini & Rodrigo Rojas                         #  
# Started: 09-Nov-2011,                                                        #
# Updates: 13-Jan-2012 ; 23-Feb-2012 ; 06-Dec-2012                             #        
################################################################################
# Purpose:                                                                     #
# This funtion read the following output files of hydroPSO:                    #
# 1) BestParameterSet.txt    : best parameter set and its corresponding        #
#                              goodness-of-fit found during the optimisation   #
# 2) Particles.txt           : parameter values and their corresponding        #
#                              goodness-of-fit for all the particles and       #
#                              iterations                                      #
# 3) Velocities.txt          : velocity values and their corresponding         #
#                              goodness-of-fit for all the particles and       #
#                              iterations                                      #
# 4) Model_out.txt           : values of the objective function / model for    #
#                              each particle and iteration                     #
# 5) ConvergenceMeasures.txt : convergence measures                            #
# 6) Particles_GofPerIter.txt: goodness-of-fit only for all the particles      #
#                              during all the iterations                       #
################################################################################
# Output:                                                                      #
# a list with the following elements                                           #
# 1) best.param: numeric with the best parameter set                           #
# 2) best.gof  : numeric with the best fitness value of the objective function #
# 3) params    : data.frame with all the parameter sets tested during the      #
#                optimisation                                                  #
# 4) gofs      : numeric with all the fitness values computed during the       #
#                optimisation (each elment in 'gofs' corresponds to one row of #
#                'params')                                                     #                                                                
# 5) model.values:numeric or matrix/data.frame with the values of the objective#
#                 function / model for each particle and iteration             #
# 6) convergence.measures: matrix/data.frame with the convergence measues.     #
#                          See'read_convergence.R' function                    #
# 7) Particles.GofPerIter: matrix/data.frame with the goodness-of-fit only for #
#                          all the particles during all the iterations         #
################################################################################

read_results <- function(drty.out="PSO.out",
                         MinMax=NULL, 
                         beh.thr=NA, 
                         modelout.cols=NULL, # 'read_out' , 'plot_out' argument
                         nsim=NULL,          # 'read_out' argument
                         verbose=TRUE) {

   ########################       Checkings      ###############################
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
   
   # Full path to 'drty.out'
   if (basename(drty.out) == drty.out) 
     drty.out <- paste(getwd(), "/", drty.out, sep="")
   
   ######################## I) Reading #########################################

   # Working directory
   wd.bak <- getwd()   
   setwd(drty.out)
   on.exit(setwd(wd.bak))

   if (verbose) message("[                                               ]")
   if (verbose) message("[         Reading output files ...              ]")
   if (verbose) message("[                                               ]")
   
   
   # 1) File "BestParameterSet.txt"       
   #best       <- read_best()     
   #best.param <- best[["best.param.values"]]  
   #best.gof   <- best[["best.param.gof"]]
   
   # 2) File "Particles.txt"
   part       <- read_particles(MinMax=MinMax, beh.thr=beh.thr, plot=FALSE)     
   params     <- part[["part.params"]]  
   gofs       <- part[["part.gofs"]] 
   best.param <- part[["best.param"]]  
   best.gof   <- part[["best.gof"]]
   
   # 3) File "Velocities.txt"        
   veloc      <- read_velocities(MinMax=MinMax, beh.thr=beh.thr, plot=FALSE)     
   velocities <- veloc[["velocities"]]  
   
   # 4) File "Model_out.txt"        
   out <- read_out(modelout.cols=modelout.cols, MinMax=MinMax, 
                   beh.thr=beh.thr, verbose=verbose, plot=FALSE)     
   model.values <- out[["model.values"]]
   model.best   <- out[["model.best"]]    
   model.obs    <- out[["model.obs"]]
   
   # 5) File "ConvergenceMeasures.txt"        
   conv <- read_convergence(MinMax=MinMax, beh.thr=beh.thr, plot=FALSE)   
   
   # 6) File "Particles_GofPerIter.txt'        
   Particles.GofPerIter <- read_GofPerParticle(plot=FALSE)
   
   # Creating the final output
   out <- list(best.param=best.param,
               best.gof=best.gof,
               params=params,
               gofs=gofs,
               velocities=velocities,
               model.values=model.values,
               model.best=model.best,
               model.obs=model.obs,
               convergence.measures=conv,
               part.GofPerIter=Particles.GofPerIter
               )
   
   # Output
   return(out)   
         
} # 'read_results' END
