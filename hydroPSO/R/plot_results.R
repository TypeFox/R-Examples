# File plot_results.R
# Part of the hydroPSO R package, http://www.rforge.net/hydroPSO/ ; 
#                                 http://cran.r-project.org/web/packages/hydroPSO
# Copyright 2011-2012 Mauricio Zambrano-Bigiarini & Rodrigo Rojas
# Distributed under GPL 2 or later

################################################################################
#                         'plot_results'                                       # 
################################################################################
# Author : Mauricio Zambrano-Bigiarini & Rodrigo Rojas                         #  
# Started: 10-Nov-2011,                                                        #
# Updates: 13-Ene-2012 ; 15-Feb-2012 ; 21-Feb-2012 ; 09-Mar-2012 ; 23-Mar-2012 #       
#          11-Jun-2012 ; 12-Nov-2012 ; 06-Dec-2012                             # 
#          21-Feb-2013                                                         #
#          09-Abr-2014                                                         #
################################################################################

plot_results <- function(drty.out="PSO.out",
                         #files=c("Particles.txt", "BestParameterSet.txt", "Model_out.txt", "ConvergenceMeasures.txt", "Velocities.txt"),                         
                         
                         ### plot.particles parameters ###
                         param.names,
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
                         cex.main=1.7,
                         cex.axis=1.3,
                         cex.lab=1.5,
                         #...,  
                         breaks="Scott",
                         freq=TRUE,
                         do.pairs=FALSE, 

                         #######################################################
                         #### For ECDFs of parameter values ('params2ecdf') ####
                         weights=NULL,                                                  
                         byrow=FALSE, 
                         leg.cex=1.2,
                         
                         #######################################################         
                         ## Parameters for the 3D dotty plots ('plot_NparOF') ###
                         dp3D.names="auto",
                         GOFcuts="auto",
                         colorRamp= colorRampPalette(c("darkred", "red", "orange", "yellow", "green", "darkgreen", "cyan")), 
                         alpha=0.65,   
                         points.cex=0.7, 
                         
                         #######################################################
                         # Parameters for GoF per Particle ('GofPerParticle')
                         ptype="one",                         
                         
                         #######################################################
                         # Parameters for BestSim vs Obs ('read_out')
                         nsim=NULL,
                         
                         #######################################################
                         # Parameters for BestSim vs Obs ('plot_out')
                         modelout.cols=NULL,
                         ftype="o", 
                         FUN=mean,                           
                         #### OPTIONS for ('plot_out') #####
                         quantiles.desired= c(0.05, 0.5, 0.95),
                         quantiles.labels= c("Q5", "Q50", "Q95"),
                         #leg.pos="bottomright",
                         
                         #######################################################
                         # parameter for convergence plot ('plot_convergence') #
                         legend.pos="topright",                            
                         
                         #######################################################
                         #####################     PNG options     #############
                         do.png=FALSE,
                         png.width=1500,
                         png.height=900,
                         png.res=90,
                         #png.drty="pngs",
                         dotty.png.fname="Params_DottyPlots.png",
                         hist.png.fname ="Params_Histograms.png",
                         bxp.png.fname="Params_Boxplots.png",
                         ecdf.png.fname ="Params_ECDFs.png",
                         pruns.png.fname="Params_ValuesPerRun.png",
                         dp3d.png.fname ="Params_dp3d.png",                                 
                         pairs.png.fname="Params_Pairs.png",
                         part.png.fname ="Particles_GofPerIter.png",
                         vruns.png.fname="Velocities_ValuePerRun.png",
                         modelout.best.png.fname="ModelOut_BestSim_vs_Obs.png",
                         modelout.quant.png.fname="ModelOut_Quantiles.png",
                         conv.png.fname ="ConvergenceMeasures.png",
                         
                         verbose=TRUE
                         ) {
   
   ######################## I) Reading #########################################
   
   # Full path to 'drty.out'
   if (basename(drty.out) == drty.out) 
     drty.out <- paste(getwd(), "/", drty.out, sep="")
     
   # Checking 'drty.out' if necessary
   if ( !file.exists(drty.out) )
     stop("Invalid argument: the directory '", drty.out, "' does not exist !")

   # PNG directory
   png.drty <- "pngs"
     
   # Full path to 'png.drty'
   if (basename(png.drty) == png.drty) 
     png.drty <- paste(drty.out, "/", png.drty, sep="")
     
   # Creating 'png.drty' if necessary
   if ( do.png & (!file.exists(png.drty)) )
     dir.create(png.drty)

   # Adding path to PNG files, if necessary
   if (basename(dotty.png.fname) == dotty.png.fname) 
      dotty.png.fname <- paste(png.drty, "/", dotty.png.fname, sep="")
   if (basename(hist.png.fname) == hist.png.fname)   
      hist.png.fname <- paste(png.drty, "/", hist.png.fname, sep="")
   if (basename(bxp.png.fname) == bxp.png.fname)   
      bxp.png.fname <- paste(png.drty, "/", bxp.png.fname, sep="")
   if (basename(ecdf.png.fname) == ecdf.png.fname)   
      ecdf.png.fname <- paste(png.drty, "/", ecdf.png.fname, sep="")
   if (basename(pruns.png.fname) == pruns.png.fname) 
      pruns.png.fname <- paste(png.drty, "/", pruns.png.fname, sep="")
   if (basename(dp3d.png.fname) == dp3d.png.fname)   
      dp3d.png.fname <- paste(png.drty, "/", dp3d.png.fname, sep="")
   if (basename(pairs.png.fname) == pairs.png.fname) 
      pairs.png.fname <- paste(png.drty, "/", pairs.png.fname, sep="")
   if (basename(part.png.fname) == part.png.fname)   
      part.png.fname <- paste(png.drty, "/", part.png.fname, sep="")
   if (basename(vruns.png.fname) == vruns.png.fname) 
      vruns.png.fname <- paste(png.drty, "/", vruns.png.fname, sep="")
   if (basename(modelout.best.png.fname) == modelout.best.png.fname) 
      modelout.best.png.fname <- paste(png.drty, "/", modelout.best.png.fname, sep="")
   if (basename(modelout.quant.png.fname) == modelout.quant.png.fname) 
      modelout.quant.png.fname <- paste(png.drty, "/", modelout.quant.png.fname, sep="")
   if (basename(conv.png.fname) == conv.png.fname)   
      conv.png.fname <- paste(png.drty, "/", conv.png.fname, sep="")

   #############################################################################
   # 1.1) Reading all the results of hydroPSO
   res <- read_results(drty.out=drty.out, MinMax=MinMax, beh.thr=beh.thr, 
                       modelout.cols=modelout.cols, verbose=verbose)
   #############################################################################
   
   # 1.2) Assignments
   best.param           <- res[["best.param"]]
   best.gof             <- res[["best.gof"]]
   params               <- res[["params"]]
   gofs                 <- res[["gofs"]]
   velocities           <- res[["velocities"]]
   model.values         <- res[["model.values"]]
   model.best           <- res[["model.best"]]
   model.obs            <- res[["model.obs"]]
   
   # If params.sub was provided
   if (!missing(param.names)) {
   
     # Number of parameters that will be analysed
     npar <- length(param.names)
    
     # Checking 'param.names'
     for ( i in 1:npar) {
       if ( !(param.names[i] %in% colnames(params)) )
         stop("Invalid argument: The field '", param.names[i], "' does not exist in 'params' !")
       } # IF end
       
     # Subsetting
     params     <- params[, param.names]
     velocities <- velocities[, param.names]
       
   } # IF end

   # Conv
   conv                 <- res[["convergence.measures"]]
   #GbestPerIter        <- res[["GbestPerIter"]]
   #GbestRate           <- res[["GbestRate"]]
   #BestFitPerIter      <- res[["BestFitPerIter"]]
   #normSwarmRadius     <- res[["normSwarmRadius"]]
   #GbestPbestRatio     <- res[["GbestPbestRatio"]]
   Particles.GofPerIter <- res[["part.GofPerIter"]]   
   
   
   ######################## II) Plotting #######################################  
   if (verbose) message("[                                               ]")
   if (verbose) message("[                  Plotting ...                 ]")
   if (verbose) message("[                                               ]")  
  
   # 2.1) Plotting parameter values: 
   #      1) Dotty Plots, 
   #      2) Histograms,
   #      3) Boxplots 
   #      4) Correlation Matrix (optional)
   #      5) Empirical CDFs
   #      6) Parameter Values Against Number of Model Evaluations 
   #      7) (pseudo)3D dotty plots
   plot_particles(#### 'plotparam' parameters ####
                  params=params, 
                  gofs= gofs,
                  gof.name=gof.name,
                  MinMax=MinMax, 
                  beh.thr=NA, 
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

                  #####################################################
                  # For ECDFs of parameter values ('params2ecdf')
                  weights=weights,                                                  
                  byrow=byrow, 
                  leg.cex=leg.cex,
                  
                  #### Parameters for the 3D dotty plots ('plot_NparOF') ####
                  dp3D.names=dp3D.names,
                  GOFcuts=GOFcuts,
                  colorRamp= colorRamp, 
                  alpha=alpha,
                  points.cex=points.cex, 
                  
                  #### parameter for convergence plot ('plot_convergence') #####
                  legend.pos=legend.pos,  
                  
                  #####################     PNG options     ####################
                  do.png=do.png,
                  png.width=png.width,
                  png.height=png.height,
                  png.res=png.res,
                  dotty.png.fname=dotty.png.fname,
                  hist.png.fname=hist.png.fname,
                  bxp.png.fname=bxp.png.fname,
                  ecdf.png.fname=ecdf.png.fname,                  
                  runs.png.fname=pruns.png.fname,
                  dp3d.png.fname=dp3d.png.fname,
                  pairs.png.fname=pairs.png.fname 
                  )
    
   # 2.2) Plotting GoF for each particle against Number of Model Evaluations
   if (!do.png) dev.new()
   plot_GofPerParticle(x=Particles.GofPerIter,
                       ptype=ptype,
                       main=main,
                       #xlab="Number of Iterations",
                       nrows="auto",
                       cex=cex,
                       cex.main=cex.main,
                       cex.axis=cex.axis,
                       cex.lab=cex.axis,      
                       col=rainbow(ncol(Particles.GofPerIter)),
                       lty=3,
                       ylim=NULL,
                       verbose=FALSE,
                       # PNG options
                       do.png=do.png,
                       png.width=png.width,
                       png.height=png.height,
                       png.res=png.res,
                       png.fname=part.png.fname  
                       )
    
    
    # 2.3) velocity values vs Number of Model Evaluations
    msg <- "[ Plotting velocity values vs Number of Model Evaluations"
    if (do.png) msg <- paste(msg, " into '", basename(vruns.png.fname), sep="")
    msg <- paste(msg, "' ...]", sep="")
    if (verbose) message(msg) 
   
    if (!do.png) dev.new()
    param.names <- colnames(params)
    plot_ParamsPerIter(params=velocities, 
                       param.names=paste("Vel.",param.names,sep=""),
                       main=main,
                       xlab="Model evaluations",
                       nrows=nrows,
                       cex=cex,
                       cex.main=cex.main,
                       cex.axis=cex.axis,
                       cex.lab=cex.lab,
                       col=rainbow(ncol(velocities)),
                       lty=3,
                       verbose=FALSE,
                       # PNG options
                       do.png=do.png,
                       png.width=png.width,
                       png.height=png.height,
                       png.res=png.res,
                       png.fname=vruns.png.fname  
                       )  

   # 2.4) Plotting Sim vs Obs
   obs.is.zoo <- FALSE
   
   if ( (length(model.best) > 1) & is.numeric(model.obs) ) {
     
     L     <- nchar(modelout.best.png.fname)
     fname <- substr(modelout.best.png.fname, 1, L-4)     
     
     if ( is.zoo(model.obs) ) { # zoo::is.zoo
       if(class(time(model.obs))=="Date") obs.is.zoo <- TRUE
     } # IF end
     
     # 2.4.1) Correlation between Best Sim and Obs
     if ( obs.is.zoo ) {         
       fname2 <- paste(fname, "-Corr.png", sep="")
     } else fname2 <- modelout.best.png.fname

     if (!do.png) dev.new() 
     plot_out(sim=model.best, 
              obs=model.obs, 
              dates=NULL, 
              ptype="corr", 
              MinMax=MinMax, 
              ftype=ftype, 
              FUN=FUN,    
              verbose=TRUE,
              
              ####
              main=main,
              
              leg.cex=leg.cex,
              #leg.pos="bottomright",              
              cex.axis=cex.axis, 
              cex.main=cex.main, 
              cex.lab=cex.lab,
              #### PNG options ### 
              do.png=do.png, 
              png.width=png.width,
              png.height=png.height, 
              png.res=png.res,
              png.fname=fname2
             )
              
     # 2.4.2) ggof between Best Sim and Obs        
     if( obs.is.zoo ) {     
        fname2 <- paste(fname, "-ggof.png", sep="")
        if (!do.png) dev.new()
        plot_out(sim=model.best, 
                 obs=model.obs, 
                 dates=NULL, 
                 ptype="ts", 
                 MinMax=MinMax, 
                 ftype=ftype, 
                 FUN=FUN,    
                 verbose=TRUE,
                 
                 ####
                 main=main,
                 
                 leg.cex=leg.cex,
                 #leg.pos="bottomright",              
                 cex.axis=cex.axis, 
                 cex.main=cex.main, 
                 cex.lab=cex.lab,
                 #### PNG options ### 
                 do.png=do.png, 
                 png.width=png.width,
                 png.height=png.height, 
                 png.res=png.res,
                 png.fname=fname2
                )
     } # IF end
     
   } # IF end
                           
   # 2.5) Plotting ECDFs for model's output OR ECDFS for quantiles of model's output
   if (!do.png) dev.new()
   if( obs.is.zoo ) { 
        plot_out(sim=model.values, 
                 obs=model.obs, 
                 dates=NULL, 
                 ptype="quant2ecdf", 
                 MinMax=MinMax, 
                 ftype=ftype, 
                 FUN=FUN,    
                 verbose=TRUE,
                 ####
                 weights=weights,
                 byrow=TRUE, 
                 quantiles.desired= quantiles.desired,
                 quantiles.labels= quantiles.labels,
                 #main=NULL,
                 ylab="Probability",
                 col="blue",
                 
                 ####
                 main=main,
                 
                 leg.cex=leg.cex,
                 #leg.pos="bottomright",
                 cex.axis=cex.axis, 
                 cex.main=cex.main, 
                 cex.lab=cex.lab,
                 #### PNG options ### 
                 do.png=do.png, 
                 png.width=png.width*0.67,
                 png.height=png.height*0.67, # For
                 png.res=png.res,
                 png.fname=modelout.quant.png.fname
                )
     } else {   
        # ecdf only
        plot_out(sim=model.values, 
                 obs=model.obs, 
                 dates=NULL, 
                 ptype="ecdf", 
                 MinMax=MinMax, 
                 ftype=ftype, 
                 FUN=FUN,    
                 verbose=TRUE,
                 ####
                 weights=weights,
                 byrow=TRUE, 
                 quantiles.desired= quantiles.desired,
                 quantiles.labels= quantiles.labels,
                 #main=NULL,
                 ylab="Probability",
                 col="blue",
                 
                 leg.cex=leg.cex,
                 #leg.pos="bottomright",
                 cex.axis=cex.axis, 
                 cex.main=cex.main, 
                 cex.lab=cex.lab,
                 #### PNG options ### 
                 do.png=do.png, 
                 png.width=png.width,
                 png.height=png.height, 
                 png.res=png.res,
                 png.fname=modelout.quant.png.fname
                )
                # IF end
         } # ELSE end
   
    
   # 2.6) Plotting Convergence Measures (Gbest and normSwarmRadius) vs Iteration Number
   if (!do.png) dev.new()
   plot_convergence(x=conv,
                    legend.pos=legend.pos,
                    verbose=TRUE,
                    #cex=1, 
                    cex.main=cex.main, 
                    cex.axis=cex.axis, 
                    cex.lab=cex.lab, 
                    # PNG options
                    do.png=do.png,
                    png.width=png.width,
                    png.height=png.height,
                    png.res=png.res,
                    png.fname=conv.png.fname,
                    ) 
    
   # 3) END               
   if (verbose) message("[                                               ]")
   if (verbose) message("[             Plots are finished !!             ]")
   if (verbose) message("[                                               ]")  
         
} # 'plot_results' END
