# File plot_particles.R
# Part of the hydroPSO R package, http://www.rforge.net/hydroPSO/ ; 
#                                 http://cran.r-project.org/web/packages/hydroPSO
# Copyright 2011-2012 Mauricio Zambrano-Bigiarini & Rodrigo Rojas
# Distributed under GPL 2 or later

################################################################################
#                          'plot_particles'                                    # 
# Author : Mauricio Zambrano-Bigiarini & Rodrigo Rojas                         #  
# Started: 08-Nov-2011,                                                        #
# Updates: 02-Feb-2012 ; 17-Feb-2012 ; 09-Mar-2012 ; 12-Nov-2012 ; 19-Nov-2012 #   
#          27-Feb-2013                                                         #   
#          09-Abr-2014                                                         #  
################################################################################
# This function plots the contents of the 'Particles.txt' ouput file of        #
# hydroPSO, with the position and fitness value of all the particles in the    #
# swarm for all the iterations.                                                #
# The following plots are produced:                                            #
# 1) Dotty Plots of Parameter Values                                           #
# 2) Histograms of Parameter Values                                            #
# 3) Boxplots of Parameter Values                                              #
# 4) Correlation matrix among Parameter Values (optional)                      #
# 5) Empirical CDFs of Parameter Values                                        #
# 6) Parameter Values Against Number of Model Evaluations                      #
# 7) Plotting pseudo-3D dotty plots of Parameter Values                        #
################################################################################
                      
plot_particles <- function(#####################################################
                           #### 'plotparam' parameters ####
                           params, #parameter values
                           gofs,
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

                           #####################################################
                           # For ECDFs of parameter values ('params2ecdf')
                           weights=NULL,                                                  
                           byrow=FALSE, 
                           leg.cex=1.5,
                                       
                           #####################################################         
                           # Parameters for the 3D dotty plots ('plot_NparOF') ##
                           dp3D.names="auto",
                           GOFcuts="auto",
                           colorRamp= colorRampPalette(c("darkred", "red", "orange", "yellow", "green", "darkgreen", "cyan")), 
                           alpha=1,
                           points.cex=0.7,  
                           legend.pos="topleft", # not used yet,but included for avoiding warnings from 'plot_results'
                           verbose=TRUE,  
                           
                           #####################################################
                           #### PNG options
                           do.png=FALSE,
                           png.width=1500,
                           png.height=900,
                           png.res=90,
                           #png.drty="pngs",
                           dotty.png.fname="Params_DottyPlots.png",
                           hist.png.fname="Params_Histograms.png",
                           bxp.png.fname="Params_Boxplots.png",
                           ecdf.png.fname="Params_ECDFs.png",
                           runs.png.fname="Params_ValuesPerRun.png",
                           dp3d.png.fname="Params_dp3d.png",
                           pairs.png.fname="Params_Pairs.png"
                           ) {
                      
  # Checking that 'params' exists
  if ( missing(params) ) stop( "Missing argument: 'params'" )
  
  # number of parameters
  nparam <- ncol(params)

  # Number of parameter sets
  n <- NROW(params)

  # Checking 'gofs'
  if (missing(gofs)) {
    stop("Missing argument: 'gofs' must be provided !!" )
  } else if (length(gofs) != n)
      stop("Invalid argument: 'length(gofs) != nrow(params)' (", length(gofs), "!=", n, ") !!" )
  
  # Parameter names
  param.names <- colnames(params)
  
  # Adding the GoF to the used data.frame
  z <- cbind(params, gofs)
  colnames(z)[nparam+1] <- gof.name
  
  ##############################################################################  
  
  # 1) Plotting Dotty Plots of Parameter Values
  if (verbose) message( "                                                     ")  
  msg <- "[ Plotting dotty plots for parameter values"
  if (do.png) msg <- paste(msg, " into '", basename(dotty.png.fname), sep="")
  msg <- paste(msg, "' ... ]", sep="")
  if (verbose) message(msg)    
   
  plot_params(params=params, 
              gofs=gofs,
              ptype="dottyplot",
              param.cols=1:nparam,
              param.names=param.names,
              #of.col=ncol(z), 
              of.name=gof.name, 
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
              verbose=FALSE,
              # PNG options
              do.png=do.png,
              png.width=png.width,
              png.height=png.height,
              png.res=png.res,
              png.fname=dotty.png.fname  
              )
            
   #############################################################################
   # 2) Plotting Histograms of Parameter Values
   msg <- "[ Plotting histograms for parameter values"
   if (do.png) msg <- paste(msg, " into '", basename(hist.png.fname), sep="")
   msg <- paste(msg, "' ... ]", sep="")
   if (verbose) message(msg)    
   if (!do.png) dev.new()
   
   plot_params(params=params, 
               gofs=gofs, 
               ptype="histogram",
               param.cols=1:nparam,
               param.names=param.names,
               #of.col=ncol(z), 
               of.name=gof.name, 
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
               verbose=FALSE,
               # PNG options
               do.png=do.png,
               png.width=png.width,
               png.height=png.height,
               png.res=png.res,
               png.fname=hist.png.fname  
               )     
               
   #############################################################################
   # 3) Plotting boxplots of Parameter Values
   msg <- "[ Plotting boxplots for parameter values"
   if (do.png) msg <- paste(msg, " into '", basename(bxp.png.fname), sep="")
   msg <- paste(msg, "' ... ]", sep="")
   if (verbose) message(msg)    
   if (!do.png) dev.new()
   
   plot_params(params=params, 
               gofs=gofs, 
               ptype="boxplot",
               param.cols=1:nparam,
               param.names=param.names,
               #of.col=ncol(z), 
               of.name=gof.name, 
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
               verbose=FALSE,
               # PNG options
               do.png=do.png,
               png.width=png.width,
               png.height=png.height,
               png.res=png.res,
               png.fname=bxp.png.fname  
               )     
            
   #############################################################################
   # 4) Plotting Correlation Matrix of Parameter Values (with hydroTSM::hydropairs)  
   if (do.pairs) 
     if ( require(hydroTSM) ) {
       msg <- "[ Plotting correlation matrix for parameter values"
       if (do.png) msg <- paste(msg, " into '", basename(pairs.png.fname), sep="")
       msg <- paste(msg, "' ... ]", sep="")
       if (verbose) message(msg)    
       if (!do.png) dev.new()
   
       plot_params(params=params, 
                   gofs=gofs, 
                   ptype="pairs",
                   param.cols=1:nparam,
                   param.names=param.names,
                   #of.col=ncol(z), 
                   of.name=gof.name, 
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
                   verbose=FALSE,
                   # PNG options
                   do.png=do.png,
                   png.width=png.width,
                   png.height=png.height,
                   png.res=png.res,
                   png.fname=pairs.png.fname  
                   )     
     } else warning("'hydroTSM' package is missing: Correlation among parameters was not plotted !")
   
   #############################################################################
   # 5) Empirical CDFs of Parameter Values
   msg <- "[ Plotting empirical CDFs for parameter values"
   if (do.png) msg <- paste(msg, " into '", basename(ecdf.png.fname), sep="")
   msg <- paste(msg, "' ... ]", sep="")
   if (verbose) message(msg)    
   if (!do.png) dev.new()
   
   params2ecdf(params=params, 
               param.names=param.names,               
               gofs= gofs,
               MinMax=MinMax, 
               beh.thr=beh.thr,                
               weights=weights,                                                  
               byrow=byrow, 
               plot=TRUE,
               nrows=nrows,  
               ylab="Probability",
               main=main,
               col=col,
               leg.cex=leg.cex,
               cex.axis=cex.axis, 
               cex.main=cex.main, 
               cex.lab=cex.lab,
               verbose=FALSE, 
               #...,
               # PNG options
               do.png=do.png,
               png.width=png.width,
               png.height=png.height,
               png.res=png.res,
               png.fname=ecdf.png.fname  
               )    
    
   #############################################################################
   # 6) Parameter Values Against Number of Model Evaluations
   msg <- "[ Plotting parameter values vs Number of Model Evaluations"
   if (do.png) msg <- paste(msg, " into '", basename(runs.png.fname), sep="")
   msg <- paste(msg, "' ... ]", sep="")
   if (verbose) message(msg)    
   if (!do.png) dev.new()
   
   plot_ParamsPerIter(params=params, 
                  param.names=param.names,
                  main=main,
                  xlab="Model Evaluations",
                  nrows=nrows,
                  cex=cex,
                  cex.main=cex.main,
                  cex.axis=cex.axis,
                  cex.lab=cex.lab,
                  col=rainbow(ncol(params)),
                  lty=3,
                  verbose=FALSE,
                  #...,
                  # PNG options
                  do.png=do.png,
                  png.width=png.width,
                  png.height=png.height,
                  png.res=png.res,
                  png.fname=runs.png.fname  
                  )
    
    ############################################################################
    # 7) Plotting 3D dotty plots of Parameter Values
    if (length(dp3D.names)==1) {
      if (dp3D.names == "auto") {
         if (nparam > 5)  nparam <- ceiling(nparam/2)
         params <- params[,1:nparam]       
         # Adding the GoF to the used data.frame
         z <- cbind(params, gofs)
         # number of parameters to be used in 3D dotty plots
         nparam <- ncol(params)
         colnames(z)[nparam+1] <- gof.name
         # Parameter names
         dp3D.names <- colnames(params)
      } # IF end
   } # IF end 
   msg <- "[ Plotting 3D dotty plots for parameter values"
   if (do.png) msg <- paste(msg, " into '", basename(dp3d.png.fname), sep="")
   msg <- paste(msg, "' ... ]", sep="")
   if (verbose) message(msg)    
    
    if (do.png) png(filename=dp3d.png.fname, width=png.width, height=png.height, res=png.res)
    else dev.new() 
     
    plot_NparOF(#params=z, 
                params=params,
                gofs=gofs,
                param.names=dp3D.names,
                MinMax=MinMax, 
                beh.thr=beh.thr,
                nrows=nrows,
                gof.name=gof.name, 
                #main=main,
                main="", # removed for avoiding confusion when 'main' is too long
                GOFcuts=GOFcuts,
                colorRamp= colorRamp,
                points.cex=points.cex,  
                alpha=alpha,
                #axis.rot=axis.rot,            
                verbose=FALSE
                ) 
    if (do.png) dev.off()
    
}  # 'plot_particles' END
