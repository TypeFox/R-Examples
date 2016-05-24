# File plot_out.R
# Part of the hydroPSO R package, http://www.rforge.net/hydroPSO/ ; 
#                                 http://cran.r-project.org/web/packages/hydroPSO
# Copyright 2012 Mauricio Zambrano-Bigiarini & Rodrigo Rojas
# Distributed under GPL 2 or later

################################################################################
#                            'plot_out'                                        # 
################################################################################
# Purpose: This function plot the values of the objective function / model for #
#          each particle and iteration                                         #
################################################################################
# Output : list with two elements:                                             #
#         -) sim: numeric with the best simulated values of the model /        #
#                 objective function among all the particles and iterations    #
#         -) obs: numeric with the observed values used during the optimisation#
################################################################################
# Author : Mauricio Zambrano-Bigiarini                                         #  
# Started: 03-Feb-2012                                                         #
# Updates: 15-Feb-2012 ; 22-Feb-2012 ; 28-Jun-2012 ; 26-Sep-2012               #   
#          21-Feb-2013                                                         #     
################################################################################

plot_out <- function(sim, obs, 
                     dates=NULL, 
                     ptype=c("corr", "ts", "ecdf", "quant2ecdf"), 
                     MinMax=NULL, 
                     #beh.thr=NA, 
                     ftype="o",  # OPTIONAL, only used when 'ptype=="ts"'.See [hydroGOF]{ggof}
                     FUN=mean,   # OPTIONAL, only used when 'ptype=="ts"'.See [hydroGOF]{ggof}
                     verbose=TRUE,
                     
                     #### OPTIONS for ptype="quant2ecdf' #####
                     weights=NULL,
                     byrow=TRUE, 
                     quantiles.desired= c(0.05, 0.5, 0.95),
                     quantiles.labels= c("Q5", "Q50", "Q95"),
                     main=NULL,
                     ylab="Probability",
                     col="blue",
                     
                     # General plotting options
                     leg.cex=1.2,
                     leg.pos="bottomright",
                     cex.axis=1.2, 
                     cex.main=1.2, 
                     cex.lab=1.2,
                       
                     #### PNG options ### 
                     do.png=FALSE, 
                     png.width=1500, 
                     png.height=900, 
                     png.res=90,
                     png.fname="ModelOut_vs_Obs.png") {                         
                         
  
  ##############################################################################
  # 1)                            Checkings                                    #
  ##############################################################################
  # Checking 'sim'
  if ( missing(sim) ) stop( "Missing argument: 'sim'" )

  # Checking 'obs'
  if ( missing(obs) ) stop( "Missing argument: 'obs'" ) 

  # Setting 'ptype' 
  ptype <- match.arg(ptype)       
  
  # number of model outputs for each parameter set
  ifelse( (is.matrix(sim) | is.data.frame(sim)), nouts <- ncol(sim), 
                                                 nouts <- length(sim))  
  
  # Checking dates' length
  if (!is.null(dates)) {
    if ( length(dates) != nouts) 
      stop("Invalid argument: 'length(dates) != length(sim)' ", length(dates), "!=", nouts, " !!")
  } #IF end 
  
#  # Checking 'MinMax'
#  if ( is.null(MinMax) ) {
#     stop("Missing argument: 'MinMax' must be provided !!")
#  } else if ( !(MinMax %in% c("min", "max")) )
#           stop("Invalid argument: 'MinMax' must be in c('min', 'max')")
  
  # Checking 'hydroGOF' pacakge when ptype=="ts"
  if (ptype=="ts") {
   if ( is.na( match("hydroGOF", installed.packages()[,"Package"] ) ) )
     stop("Invalid argument: You don't have the 'hydroGOF' package => You can not use 'ptype='ts' !!")
  } # IF end
           
  # Checking 'class(sim)'    
  if ( (ptype=="corr") | (ptype=="ts") ) {
    if ( !(class(sim) %in% c("numeric", "integer") ) )
      stop("Invalid argument: 'class(sim)' must be in c('numeric', 'integer') for 'ptype' in c('corr', 'ts') !!")
  } else if ( (ptype=="ecdf") | (ptype=="quant2ecdf") ) {
    # if ( (class(sim) != "matrix") & (class(sim) != "data.frame") )
    #  stop("Invalid argument: 'class(sim)' must be in c('matrix', 'data.frame') for 'ptype' in c('ecdf', 'quant2ecdf') !!")
     } # ELSe end
  
  
  ##############################################################################
  # 2)                             Plotting                                    #
  ##############################################################################

  # 1) Plotting Observed against Simulated values
  if ( (ptype=="corr") | (ptype=="ts") ) {
    msg <- "[ Plotting best simulated values vs observations"
  } else if ( (ptype=="ecdf") | (ptype=="quant2ecdf") ) {
      msg <- "[ Plotting ECDFs of simulated quantiles vs observations"
    } # ELSE end
  if (do.png) msg <- paste(msg, " into '", basename(png.fname), sep="")
  msg <- paste(msg, "' ... ]", sep="")
  if (verbose) message(msg)       

  # PNG output ?
  if (do.png) {
    if (ptype=="corr") { png(filename= png.fname, width= min(png.width, png.height), height= min(png.width, png.height), res=png.res)  
    } else if (ptype!="ecdf") png(filename= png.fname, width = png.width, height = png.height, res=png.res)  
  } #IF end
   
  if (ptype=="corr") {
  
    lm.os <-  lm(sim ~ obs)
    axis.lims <- range(sim, obs, na.rm=TRUE)
    plot(obs, sim, xlab="Observations", ylab="Best Simulation", 
         ylim=axis.lims, xlim=axis.lims, cex.axis=1.2, cex.main=1.2, cex.lab=1.2, panel.first=grid() )   
    abline(lm.os, col="grey", lty=3)
    r2 <- (cor(sim, obs, use="complete.obs"))^2
    legend("right", legend=paste("R2=", round(r2,4), sep=""),  inset=0.02 )
    
  } else if (ptype=="ts") {
      if (is.null(main)) main="Observed vs 'best' Simulation"
       
      if (!is.null(dates)) {
        if (!is.zoo(obs)) obs <- zoo(obs, dates) # zoo::is.zoo ; zoo::zoo
        if (!is.zoo(sim)) sim <- zoo(sim, dates) # zoo::is.zoo ; zoo::zoo 
      } # IF end
      
      # Plotting Sim vs Obs
      library(hydroGOF)
      ggof(sim=sim, obs=obs, main=main, ftype=ftype, FUN=FUN, 
           cex.main=cex.main, cex.axis=cex.axis, cex.lab=cex.lab, leg.cex=leg.cex) # dates.fmt="%Y-%m-%d")
  
      # Computing several GoFs between the best behavioural simulation and the observed values
      #if (verbose) message("                                                  ")
      #if (verbose) message("GoFs between 'obs' and 'best' simulation: ")  
      #if (verbose) print( gof(sim=sim, obs=obs, ...) )
      #if (verbose) print( gof(sim=sim, obs=obs) )
      
    } else if (ptype=="ecdf") {
        
       params2ecdf(params=sim, 
                   #param.names=NULL,
                   weights=weights,                                                  
                   byrow=byrow, 
                   plot=TRUE,
                   obs=obs,
                   main=main,
                   nrows="auto",  
                   ylab=ylab,
                   col=col,
                   leg.cex=leg.cex,
                   leg.pos=leg.pos,
                   cex.axis=cex.axis, 
                   cex.main=cex.main, 
                   cex.lab=cex.lab,
                   verbose=TRUE, 
                   do.png=do.png,
                   png.width=png.width, 
                   png.height=png.height, 
                   png.res=png.res,
                   png.fname=png.fname)       

      } else if (ptype=="quant2ecdf") {
        
        quant2ecdf(sim=sim, 
                   weights=weights,
                   byrow=byrow, 
                   quantiles.desired= quantiles.desired,
                   plot=TRUE,
                   obs=obs, 
                   quantiles.labels= quantiles.labels,
                   main=main,
                   ylab=ylab,
                   col=col,
                   leg.cex=leg.cex,
                   leg.pos=leg.pos,
                   cex.axis=cex.axis, 
                   cex.main=cex.main, 
                   cex.lab=cex.lab,
                   verbose=verbose)
       

      } # ELSE end 
    
    if (do.png & (ptype!="ecdf")) dev.off()    
  
}  # 'plot_out' END
