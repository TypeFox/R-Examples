# File read_out.R
# Part of the hydroPSO R package, http://www.rforge.net/hydroPSO/ ; 
#                                 http://cran.r-project.org/web/packages/hydroPSO
# Copyright 2011-2012 Mauricio Zambrano-Bigiarini & Rodrigo Rojas
# Distributed under GPL 2 or later

################################################################################
#                            'read_out'                                        # 
################################################################################
# Purpose: This function reads the values of the objective function / model for#
#          each particle and iteration                                         #
################################################################################
# Output : list with two elements:                                             #
#          -) model.values: numeric or matrix/data.frame with the values of the#
#                           objective function / model for each particle and   #
#                           iteration                                          #
#          -) model.gofs  : numeric vector with the goodness-of-fit value for  #
#                           each particle and iteration                        #
#          -) model.best  : numeric with the best model output                 #
#          -) model.obs   : numeric with the observed values used for hydroPSO,#
#                           or NA if they were not provided                    #
################################################################################
# Author : Mauricio Zambrano-Bigiarini & Rodrigo Rojas                         #  
# Started: 08-Nov-2011,                                                        #
# Updates: 27-Jan-2012 ; 02-Feb-2012 ; 15-Feb-2012 ; 15-Feb-2012 ; 22-Feb-2012 #  
#          23-Mar-2012 ; 06-Dec-2012                                           # 
#          29-May-2013                                                         #  
#          09-Abr-2014                                                         #  
################################################################################
# Columns in 'of_out' are:
# Iter         : integer, with the iteration number for each row of the file
# Part         : integer, with the particle number for each row of the file
# GoF          : numeric, with the goodness of fit value for each particle and iteration
# Model_Output : numeric vector (one or more values for each row), with the model output for each particle and iteration

read_out <- function(file="Model_out.txt", 
                     modelout.cols=NULL, 
                     nsim=NULL,
                     obs, 
                     MinMax=NULL, 
                     beh.thr=NA, 
                     verbose=TRUE,
                     #### Plotting arguments ####
                     plot=TRUE,
                     ptype=c("corr", "ts", "ecdf", "quant2ecdf"), 
                     ftype="dm", # OPTIONAL, only used when 'ptype=="ts"'.See [hydroGOF]{ggof}
                     FUN=mean,   # OPTIONAL, only used when 'ptype=="ts"'.See [hydroGOF]{ggof}                     
                     #### OPTIONS for ptype="quant2ecdf' #####
                     weights=NULL,
                     byrow=TRUE, 
                     quantiles.desired= c(0.05, 0.5, 0.95),
                     quantiles.labels= c("Q5", "Q50", "Q95"),
                     main=NULL,
                     ylab="Probability",
                     col="blue",
                     leg.cex=1.2,
                     leg.pos="bottomright",
                     cex.axis=1.2, 
                     cex.main=1.2, 
                     cex.lab=1.2,                     
                     #### PNG options ### 
                     do.png=FALSE, png.width=1500, png.height=900, png.res=90,
                     png.fname="ModelOut_vs_Obs.png" ) {                         
                         
  
  ##############################################################################
  # 1)                            Checkings                                    #
  ##############################################################################
  # Checking that 'file' exists
  if ( !file.exists(file) )
     stop( "Invalid argument value: The file '", basename(file), "' doesn't exist !!")

  # Checking 'beh.thr'
  if ( !is.na(beh.thr) ) {
     if ( is.null(MinMax) ) {
        warning("Missing argument: 'MinMax' has to be provided before using 'beh.thr' !!")  
     } # IF end
  } # IF end
  
  # Checking 'MinMax'
  if ( !is.null(MinMax) ) {
     if ( !(MinMax %in% c("min", "max")) )
        stop("Invalid argument: 'MinMax' must be in c('min', 'max')")
  } # IF end

  # Checking 'plot' and 'MinMax'
  valid.types <- c("corr", "ts") 
  if (plot &  (length(which(!is.na(match(ptype, valid.types )))) <= 0) ) {
    if (is.null(MinMax)) {
      warning("Missing argument: 'MinMax' has to be provided before plotting when 'ptype' is in c('corr', 'ts') => 'plot=FALSE'")
      plot <- FALSE
    } # IF end
  } # IF end

  ##############################################################################
  # 2)                            Reading                                      #
  ##############################################################################

  # Reading the file
  if (verbose) message( "                                                     ")  
  if (verbose) message( "[ Reading the file '", basename(file), "' ... ]" )  
  if (!missing(nsim)) {
    cnames <- paste("sim", 1:nsim, sep="") 
    cnames <- c("Iter", "Part", "GoF", cnames)      
    sim  <- read.table(file=file, header=FALSE, skip=1, fill=TRUE,  col.names=cnames)
  } else sim  <- read.table(file=file, header=FALSE, skip=1, fill=TRUE)

  # Amount of total model outputs / parameter sets 
  nouts <- nrow(sim)
  if (verbose) message( "[ Total number of parameter sets: ", nouts, " ]" )       
  
  # computing the number of columns in 'file'
  ncols <- ncol(sim)  
  
  # Getting the GoF for each particle and iteration
  gofs <- sim[, 3]

  # Getting the Model Outputs for each particle and iteration
  outputs <- sim[, 4:ncols]
  
  # Computing the number of values in each model output
  lnsim  <- NCOL(outputs)
  if (verbose) message( "[ Number of model outputs for each parameter set: ", lnsim, " ]" )    
  
  # giving more meaningful names to the model outputs
  if (lnsim > 1)
    colnames(outputs) <- paste("Sim", 1:lnsim, sep="")
  
  # If the user only wants some columns of the model output file
  if (!is.null(modelout.cols)) {
    outputs <- outputs[, modelout.cols]
    
    lnsim  <- NCOL(outputs)
    if (verbose) message( "[ Number of extracted model outputs for each parameter set: ", lnsim, " ]" ) 
  } # IF end   
  
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
    if ( is.matrix(outputs) | is.data.frame(outputs) ) {
      outputs <- outputs[beh.row.index, ]
    } else outputs <- outputs[beh.row.index]
    gofs <- gofs[beh.row.index]
   
    # Amount of behavioural parameter sets 
    nbeh <- nrow(outputs)
    if (verbose) message( "[ Number of behavioural model outputs           : ", nbeh, " ]" ) 
    
  } # IF end
  
  # Selecting best model output
  if ( !is.null(MinMax) ) {
    ifelse(MinMax=="min", best.row.index <- which.min(gofs), 
                          best.row.index <- which.max(gofs) )
    if ( is.matrix(outputs) | is.data.frame(outputs) ) {
      best <- as.numeric(outputs[best.row.index, ]) 
    } else best <- as.numeric(outputs[best.row.index])
    
    # If the user only wants some columns of the model output file
    if (!is.null(modelout.cols)) best <- best[modelout.cols]
    
  } else {
      best <- NA
      plot <- FALSE
      if (verbose) warning("[ Note: You didn't provide 'MinMax' => model.best=NA & some plots will be skipped !! ]")
    } # ELSE end

  ##############################################################################
  # 3)                        Getting Observed Values                          #
  ##############################################################################

  if (missing(obs)) {
    fname <- "Observations.txt"
    if (file.exists(fname)) {
      if ( !is.na( match("zoo", installed.packages()[,"Package"] ) ) ) {
        obs   <- read.zoo(fname) # zoo::read.zoo
        dates <- time(obs)
        # If the observed data do not have dates, they are transformed into numeric
        if (class(dates)!="Date") obs <- coredata(obs)
      } else {
        obs <- read.table(fname, header=FALSE, skip=0)  
        # skippping the column with the index
        obs <- obs[,1]
      } # ELSE end 
      
      # If the user only wants some columns of the model output file
      if (!is.null(modelout.cols)) {
        obs   <- obs[modelout.cols]
        dates <- dates[modelout.cols] 
      } # IF end
  
    } else obs <- NA
  } # IF end

  # Checking length(obs)
  if ( !is.null(obs) & is.numeric(obs) ) {
    if ( length(obs) != lnsim ) 
      stop("Invalid argument: 'length(obs) != ncol(sims)' ", length(obs), "!=", lnsim, " !!")
  } # IF end
  
  ##############################################################################
  # 4)                            Output                                       #
  ##############################################################################  
  # creating the output
  out <- list(model.values=outputs, model.gofs=gofs, model.best=best, model.obs=obs)

  ##############################################################################
  # 5)                            Plotting                                     #
  ##############################################################################
  if ( plot & is.numeric(obs) ) {
  
       plot_out(sim=best, 
                obs=obs, 
                dates=dates, 
                ptype=ptype, 
                MinMax=MinMax, 
                ftype=ftype, # OPTIONAL, only used when 'ptype=="ts"'.See [hydroGOF]{ggof}
                FUN=FUN,     # OPTIONAL, only used when 'ptype=="ts"'.See [hydroGOF]{ggof}
                verbose=verbose,
                #### OPTIONS for ptype="quant2ecdf' #####
                weights=weights,
                byrow=byrow, 
                quantiles.desired= quantiles.desired,
                quantiles.labels= quantiles.labels,
                main=main,
                ylab=ylab,
                col=col,
                leg.cex=leg.cex,
                leg.pos=leg.pos,
                cex.axis=cex.axis, 
                cex.main=cex.main, 
                cex.lab=cex.lab,                     
                #### PNG options ### 
                do.png=do.png, 
                png.width=png.width, 
                png.height=png.height, 
                png.res=png.res,
                png.fname=png.fname)
  } # IF end
  
  return(out)
  
}  # 'read_out' END
