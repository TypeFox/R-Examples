# File hydromod.R
# Part of the hydroPSO R package, http://www.rforge.net/hydroPSO/ ; 
#                                 http://cran.r-project.org/web/packages/hydroPSO
# Copyright 2010-2013 Mauricio Zambrano-Bigiarini & Rodrigo Rojas
# Distributed under GPL 2 or later

################################################################################
#                           hydromod                                           #
################################################################################
# Purpose    : To run a hydrological/environmental model, and get a            #
#              goodness-of-fit value by comparing the simulated values against # 
#              observations                                                    #
################################################################################
# Output     : A list of two elements:                                         #
#              1) sim: simulated values obtained by running the hydrological   #
#                      model                                                   #
#              2) GoF: goodness-of fit of the simualted values against observed#
#                      ones, by using THE USER-DEFINED 'gof' measure           # 
################################################################################
# Author : Mauricio Zambrano-Bigiarini                                         #
# Started: 14-Dec-2010 at JRC Ispra                                            #
# Updates: 20-Dec-2010                                                         #
#          19-Jan-2011 ; 22-Jan-2011 ; 02-Feb-2011 ; 11-May-2011               #
#          13-Jan-2012 ; 16-Jan-2012 ; 23-Jan-2012 ; 02-May-2012 ; 12-Oct-2012 #
#          15-Oct-2012 ; 22-Nov-2012                                           #
#          29-May-2013                                                         #
################################################################################
hydromod <- function(
                     param.values,                 # Numeric vector with the paramter values that will be used in the input files of the hydrological model
                     param.files="ParamFiles.txt", # Character, with the name of the file (with full path) that stores the name of the files that have to be modified for each parameter 
                     model.drty=getwd(),           # Character, with the path of the directory that stores the exe file of the hydrological model and ALL the input files required for the simulation
                     exe.fname,
                     stdout=FALSE,                 # a logical (not NA) indicating whether messages written to ‘stdout’ should be sent or not. See '?system2'
                     stderr="",                    # a logical (not NA) indicating whether messages written to ‘stderr’ should be sent or not. See '?system2'
                     verbose= FALSE,               # logical, indicating if progress messages have to be printed during the simulations. If \code{verbose=TRUE}, the following messages will appear: i)parameter values for each particle; (ii) model execution; iii) extraction of simulated values; and iv) computation of the goodness-of-fit measure.
                                          
                     ###############################
                     ###############################
                     # function for reading the ouput
                     # of the hydrological model   #
                     ###############################
                     out.FUN,                      # Character, with the name of a valid R function to be used for reading the outputs of te model in 'out.fname' and transforming them into a zoo object
                     out.FUN.args,                 # list, with the arguments to be passed to 'out.FUN'
                     
                     ###############################
                     ###############################
                     # function for comparing simulated
                     # values with observations    #
                     ###############################
                     gof.FUN,             # name of the GoF function used to compare simulations and observations                     
                     gof.FUN.args=list(), # list, with additional arguments to be passed to 'out.FUN' (additional to 'sim' and 'obs')
                     gof.Ini,             # OPTIONAL. Starting date for the the goodness-of-fit function to be used for comparing observed against simulated values.
                                          # It is used to subset 'obs' (if necessary), AND to define the time period  to compare simulated ith observed values
                     gof.Fin,             # OPTIONAL. Ending date for the the goodness-of-fit function to be used for comparing observed against simulated values.
                                          # It is used to subset 'obs' (if necessary), AND to define the time period  to compare simulated ith observed values
                     date.fmt="%Y-%m-%d", # Format used to define 'gof.Ini", "gof.Fin"
                     obs,                 # zoo object with the observed values  
                     
                     ##############################
                     # OPTIONAL. for plotting     #
                     # sim vs obs                 #
                     do.png=FALSE, # logical, indicating if a png image with the results of the 'ggof' function has to be produced
                     png.fname,    # OPTIONAL, only used when \code{do.png=TRUE}. Name of the PNG file to be produced. By default its name is 'Obs_vs_Sim.png', within the 'model.drty'
                     width = 1200, # OPTIONAL, only used when \code{do.png=TRUE}. Widht of the output PNG image
                     height = 600, # OPTIONAL, only used when \code{do.png=TRUE}. Height of the output PNG image
                     res=90,       # OPTIONAL, only used when \code{do.png=TRUE}. Resolution of the output PNG image
                     
                     main,         # OPTIONAL, only used when \code{do.png=TRUE}. character, representing the main title of the plot comparing observed and simulated values
                     leg.cex=1.2, 
                     tick.tstep= "auto", # character, indicating the time step that have to be used for putting the ticks on the time axis. See 'ggof' 
                     lab.tstep= "auto",  #character indicating the time step that have to be used for putting the labels on the time axis. See 'ggof'  
                     lab.fmt=NULL        # Character, with the format to be used for drawing the time axis.                         
                                                   
                     ###############################                     
                     #...                # Addtional arguments to be passed to                    
                     
                     ) {
                     
  ##############################################################################
  # 0)                             Checkings                                   #
  ##############################################################################
    
  # Verifying 'param.values'
  #if (missing(param.values))
  #  stop( "Missing argument: 'param.values' has to be given !")
  
  if (!file.exists(param.files))
    stop( "Invalid argument: the file '", param.files, "' doesn't exist!" )
  
  if ( missing(exe.fname) )
    stop( "Missing argument: 'exe.fname'" )
    
  if ( missing(out.FUN) ) {
    stop( "Missing argument: 'out.FUN'" )
  } else out.FUN <- match.fun(out.FUN)
  
  if ( missing(gof.FUN) ) {
    stop( "Missing argument: 'gof.FUN'" )
  } else {gof.name <- gof.FUN
          gof.FUN <- match.fun(gof.FUN)
         } # ELSE end

  setwd(model.drty)
  
  if (!file.exists(exe.fname))
    stop( "Invalid argument: the file '", exe.fname, "' does not exist!" )

  if ( sessionInfo()[[1]]$os != "linux-gnu") {
    dot.pos   <- which(strsplit(exe.fname, split=character(0))[[1]] == ".")
    ext       <- substr(exe.fname,start=dot.pos+1, stop=nchar(exe.fname)) 
    if (ext=="exe") exe.fname <- substr(exe.fname,start=1, stop=dot.pos-1) 
  } # IF end

  if (missing(obs)) stop("Missing argument: 'obs' must be provided")
    
  ##############################################################################  
  # 1)     New parameter values -> input files of the hydrological model       #
  ##############################################################################  
  if (verbose) message("                                           ")
  if (verbose) message("===========================================")
  if (verbose) message("[ 1) Writing new parameter values ...     ]")
  if (verbose) message("===========================================")
  ParameterValues2InputFiles(NewValues=param.values, ParamFiles.fname=param.files,
                             verbose=verbose)                                       

  ##############################################################################
  # 2)                       Running the hydrological model                    #
  ##############################################################################
  
  if (verbose) message("===========================================")
  if (verbose) message("[ 2) Running the model ...                ]")
  if (verbose) message("===========================================")
  system2(exe.fname, stdout=stdout, stderr=stderr)
  
  ##############################################################################
  # 3)                     Extracting simulated values                         #                                 
  ##############################################################################
  
  if (verbose) message("===========================================")
  if (verbose) message("[ 3) Extracting simulated values ...      ]")
  if (verbose) message("===========================================")

  out.FUN.argsDefaults <- formals(out.FUN)
  out.FUN.args         <- modifyList(out.FUN.argsDefaults, out.FUN.args) 
  sim                  <- do.call(out.FUN, as.list(out.FUN.args))   

  ##############################################################################
  # 4)                     Goodness of fit                                     #                                 
  ##############################################################################
  
  if (verbose) message("===========================================")
  if (verbose) message("[ 4) Computing the goodness of fit ...    ]")
  if (verbose) message("===========================================")
  
  # In case different dates are given for 'sim', 'obs', 'gof.Ini', gof.Fin', 
  # 'sim' and 'obs' are subset to the time period selected for the GoF function 
  
  if ( !missing(gof.Ini) | !missing(gof.Fin) ) 
    ifelse ( grepl("%H", date.fmt, fixed=TRUE) | grepl("%M", date.fmt, fixed=TRUE) |
             grepl("%S", date.fmt, fixed=TRUE) | grepl("%I", date.fmt, fixed=TRUE) |
             grepl("%p", date.fmt, fixed=TRUE) | grepl("%X", date.fmt, fixed=TRUE),
             subdaily <- TRUE, subdaily <- FALSE )
  
  if (!missing(gof.Ini) ) {
      ifelse(subdaily, gof.Ini <- as.POSIXct(gof.Ini, format=date.fmt),
                       gof.Ini <- as.Date(gof.Ini, format=date.fmt) )
                   
      if (!is.zoo(obs)) { # zoo::is.zoo
        stop( "Invalid argument: 'obs' must be a zoo or xts object to use 'gof.Ini' !")
      } else obs <- window(obs, start=gof.Ini) 
      if (!is.zoo(sim)) { # zoo::is.zoo
        stop( "Invalid argument: 'sim' must be a zoo or xts object to use 'gof.Ini' !")
      } else sim <- window(sim, start=gof.Ini)
                    
  } # IF end 
   
  if (!missing(gof.Fin) ) {
    ifelse(subdaily, gof.Fin <- as.POSIXct(gof.Fin, format=date.fmt),
                     gof.Fin <- as.Date(gof.Fin, format=date.fmt) )
                       
    if (gof.Fin < gof.Ini) {
      stop( "Invalid argument: 'gof.Fin < gof.Ini' (", gof.Fin, " < ", gof.Ini, ")" )
    } else { 
             if (!is.zoo(obs)) { # zoo::is.zoo
               stop( "Invalid argument: 'obs' must be a zoo or xts object to use 'gof.Fin' !" )
             } else obs <- window(obs, end=gof.Fin)
             if (!is.zoo(sim)) { # zoo::is.zoo
               stop( "Invalid argument: 'sim' must be a zoo or xts object to use 'gof.Fin' !" )
             } else sim <- window(sim, end=gof.Fin)
           } # IF end
  } # IF end  
  
  nobs <- length(obs)
  nsim <- length(sim)
  
  if (nobs != nsim) 
    stop( "Invalid argument: number of observations != number of simulations (", nobs, "!=", nsim, ")'" )      
  
  gof.FUN.argsDefaults <- formals(gof.FUN)
  gof.FUN.args         <- modifyList(gof.FUN.argsDefaults, gof.FUN.args) 
  gof.FUN.args         <- modifyList(gof.FUN.args, list(sim=sim, obs=obs))  
  gof.value            <- do.call(gof.FUN, as.list(gof.FUN.args))   
  
  if (verbose) message("[", gof.name, "= ", round(gof.value,3), "]")
  
  ##############################################################################
  # 5)                    Creating the output object                           #                                 
  ##############################################################################
  out <- list(2)
  out[[1]] <- sim
  out[[2]] <- gof.value
  names(out) <- c("sim", "GoF")   
  
  ##############################################################################
  # 6)                    OPTIONAL: GoF plot                                    #                                 
  ##############################################################################
  
  if (do.png)  {
  
    if (verbose) message("===========================================")
    if (verbose) message("[ 5) Creating a PNG with 'sim' vs 'obs'...]")
    if (verbose) message("===========================================")
    
    if (missing(main)) main <- "Simulations vs Observations"
      
    if (missing(png.fname))
      png.fname <- paste(file.path(model.drty), "/Obs_vs_Sim.png", sep="")
    
    png(filename=png.fname, width= width, height= height, res=res)

    if ( !is.na( match("hydroGOF", installed.packages()[,"Package"] ) ) ) {
      library(hydroGOF)
      ggof(sim, obs, main=main, cex.main=1.5, leg.cex=leg.cex, tick.tstep=tick.tstep, lab.tstep=lab.tstep, lab.fmt=lab.fmt)
    } else plot_out(sim=sim, obs=obs, ptype="corr")

    dev.off()
    
    if (verbose) message("===========================================")
    if (verbose) message("[ See the file '", basename(png.fname), "' ")
    if (verbose) message("===========================================")
    
  } # IF end

  return(out)

} # 'hydromod' END
