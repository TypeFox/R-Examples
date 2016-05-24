## File verification.R
## Part of the hydroPSO R package, http://www.rforge.net/hydroPSO/ ; 
##                                 http://cran.r-project.org/web/packages/hydroPSO
## Copyright 2011-2014 Mauricio Zambrano-Bigiarini & Rodrigo Rojas
## Distributed under GPL 2 or later

################################################################################
#                           verification                                       #
################################################################################
# Purpose    : Run a hydrological model with different parameter sets (obtained#
#              during calibration) and to obtain the simulated time series for #
#              each one of the parameter sets                                  #
################################################################################
# Output     : A list of two elements:                                         #
#              ##1) sim: simulated values obtained by running the hydrological #
#                      model                                                   #
#              2) gofs    : goodness-of fit of the simulated values against    #
#                           observed ones, by using THE USER-DEFINED 'gof'     #
#                           measure                                            #
#              3) best.gof: goodness-of fit of the "best" parameter set        #
#              4) best.par: parameter values of the "best" paraemter set       # 
################################################################################
# Author  : Mauricio Zambrano-Bigiarini                                        #
# Started : 18-Jan-2011 at JRC Ispra                                           #
# Updates : 12-May-2011 ; 13-Feb-2012  ; 23-Feb-2012                           #
#           09-Abr-2014                                                        #
################################################################################
verification <- function(
                         fn="hydromod",  
                         par,
                         control=list(),
                         model.FUN=NULL,
                         model.FUN.args=list()                                
                        ) {
                     
  
  ##############################################################################
  #                            INPUT CHECKING                                  #
  ##############################################################################
        
  # Checking the name of the objective function
  if (missing(fn)) {
     stop("Missing argument: 'fn' must be provided")
  } else {
      fn.name <- fn
      fn      <- match.fun(fn)
    } # ELSE end       
        
  # Checking 'par'
  if (missing(par)) stop("Missing argument: 'par' must be provided")
        
          
  ########################################################################        
  con <- list(
                
          drty.in=getwd(),
          drty.out="verification", # Character, with the name of the directory that will store the results of the LH-OAT. 
          digits=7,
                
          gof.name="GoF",          # Character, only used for identifying the goodness-of-fit of each model run
          MinMax="max",            # Character, indicating if PSO have to find a minimum or a maximum for the objective function. \cr
                                   # Valid values are in: \code{c('min', 'max')} \cr
          do.plots=FALSE,
          write2disk=TRUE,
          verbose= TRUE            # logical, indicating if progress messages have to be printed
             )
              
  nmsC <- names(con)
 
  con[(namc <- names(control))] <- control
  
  if (length(noNms <- namc[!namc %in% nmsC])) 
    warning("[Unknown names in control: ", paste(noNms, collapse = ", "), " (not used) !]")
    
  drty.in        <- con[["drty.in"]]
  drty.out       <- con[["drty.out"]]
  digits         <- con[["digits"]]

  gof.name       <- con[["gof.name"]]
  MinMax         <- con[["MinMax"]]
  do.plots       <- as.logical(con[["do.plots"]])  
  write2disk     <- as.logical(con[["write2disk"]])
  verbose        <- as.logical(con[["verbose"]])
        
  ########################################################################
  ##################### Dummy checkings ##################################
     
  # 'hydromod' checkings
  if (fn.name=="hydromod") {

    # checking that 'model.FUN' is a valid function          
    if ( is.null(model.FUN) ) {
      stop( "'model.FUN' has to be defined !" )
    } else  {
        model.FUN.name <- model.FUN
        model.FUN      <- match.fun(model.FUN)
      } # ELSE end
  
    # checking 'out.FUN.args'
    if ( length(model.FUN.args)==0 ) {
      warning( "'model.FUN.args' is an empty list. Are you sure your model doesn't have any argument(s) ?" )
    } else {
        # Modifying the arguments of the hydrological model
        model.FUN.argsDefaults <- formals(model.FUN)
        model.FUN.args         <- modifyList(model.FUN.argsDefaults, model.FUN.args) 
      } # ELSe end
             
  } # IF end 
          
  if ( is.matrix(par) | is.data.frame(par) ) {
    # Computing 'P', the Dimension of the Solution Space
    P <- ncol(par)

    # Computing the number of parameter sets
    nparamsets <- nrow(par)

    # Meaningful name of each one of the parameters
    if (is.null(colnames(par))) {
      Parameter.names <- paste("Param", 1:nparamsets, sep="")
    } else Parameter.names <- colnames(par)  
  
  } else if (is.numeric(par)) {
      # Computing 'P', the Dimension of the Solution Space
      P <- length(par)
 
      # Computing the number of parameter sets
      nparamsets <- 1

      # Meaningful name of each one of the parameters
      if (is.null(names(par))) {
        Parameter.names <- paste("Param", 1:P, sep="")
      } else Parameter.names <- names(par)    

    } else stop("Invalid argument: 'class(par)' must be in c('numeric', 'matrix', 'data.frame')")
   
  if (verbose) message( "[ Number of parameter set read: ", nparamsets, " ]" )   
  if (verbose) message( "[ Number of parameters read   : ", nparamsets, " ]" )  
  if (verbose) message( "[ Parameter names             : ", paste(Parameter.names, collapse = ", "), " ]")
  
  # If the user only provided a single parameter set, it is transofrmed into a matrix
  if (nparamsets==1) par <- matrix(par, nrow=1)
        
  # Adding the parent path of 'drty.out', if it doesn't have it
  if (drty.out == basename(drty.out) )
    drty.out <- paste( getwd(), "/", drty.out, sep="")
        
  # Verifying that 'drty.out' directory exists. IF not, it is created
  if (!file.exists(file.path(drty.out))) {
    if (write2disk) {
      dir.create(file.path(drty.out))
      if (verbose) message("                                            ")
      if (verbose) message("[ Output directory '", basename(drty.out), "' was created on: '", dirname(drty.out), "' ]") 
      if (verbose) message("                                            ")
    } # IF end
  } # IF end  
  
  
  ##############################################################################
  #                            Writing Info File
  ##############################################################################  
     
  # File 'Verification-logfile.txt' #        
  InfoTXT.fname <- paste(file.path(drty.out), "/", "Verification-logfile.txt", sep="")
  InfoTXT.TextFile  <- file(InfoTXT.fname , "w+")
  #c(isOpen(Tfile, "r"), isOpen(Tfile, "w")) # both TRUE
  writeLines("================================================================================", InfoTXT.TextFile) # writing a separation line with a carriage return
    writeLines(c("Platform             :", sessionInfo()[[1]]$platform), InfoTXT.TextFile, sep="  ")
    writeLines("", InfoTXT.TextFile) # writing a blank line with a carriage return
    writeLines(c("R version            :", sessionInfo()[[1]]$version.string), InfoTXT.TextFile, sep="  ")
    writeLines("", InfoTXT.TextFile) # writing a blank line with a carriage return
    writeLines(c("hydroPSO version     :", sessionInfo()$otherPkgs$hydroPSO$Version), InfoTXT.TextFile, sep="  ")
    writeLines("", InfoTXT.TextFile) #writing a blank line with a carriage return
    writeLines(c("hydroPSO Built       :", sessionInfo()$otherPkgs$hydroPSO$Built), InfoTXT.TextFile, sep="  ")
    writeLines("", InfoTXT.TextFile) # writing a blank line with a carriage return
    writeLines("================================================================================", InfoTXT.TextFile) # writing a separation line with a carriage return
    Time.Ini <- Sys.time()
    writeLines(c("Starting Time        :", date()), InfoTXT.TextFile, sep="  ")
    writeLines("", InfoTXT.TextFile) # writing a blank line with a carriage return
    writeLines("================================================================================", InfoTXT.TextFile) # writing a separation line with a carriage return
    writeLines(c("drty.in              :", drty.in), InfoTXT.TextFile, sep=" ") 
    writeLines("", InfoTXT.TextFile) # writing a blank line with a carriage return
    writeLines(c("drty.out             :", drty.out), InfoTXT.TextFile, sep=" ") 
    writeLines("", InfoTXT.TextFile) # writing a blank line with a carriage return
    writeLines(c("digits               :", digits), InfoTXT.TextFile, sep=" ") 
    writeLines("", InfoTXT.TextFile) # writing a blank line with a carriage return
    writeLines(c("gof.name             :", gof.name), InfoTXT.TextFile, sep=" ") 
    writeLines("", InfoTXT.TextFile) # writing a blank line with a carriage return
    writeLines(c("do.plots             :", do.plots), InfoTXT.TextFile, sep=" ") 
    writeLines("", InfoTXT.TextFile) # writing a blank line with a carriage return
    writeLines(c("write2disk           :", write2disk), InfoTXT.TextFile, sep=" ") 
    writeLines("", InfoTXT.TextFile) # writing a blank line with a carriage return
    writeLines(c("verbose              :", verbose), InfoTXT.TextFile, sep=" ") 
    writeLines("", InfoTXT.TextFile) # writing a blank line with a carriage return
    if (fn.name=="hydromod") {
      try(writeLines(c("hydromod function    :", model.FUN.name), InfoTXT.TextFile, sep=" ") , TRUE)
      writeLines("", InfoTXT.TextFile) # writing a blank line with a carriage return
      writeLines(c("hydromod args        :"), InfoTXT.TextFile, sep=" ") 
      writeLines("", InfoTXT.TextFile) # writing a blank line with a carriage return
      for ( i in 1:length(model.FUN.args) ) {              
        arg.name  <- names(model.FUN.args)[i]
        arg.name  <- format(paste("  ", arg.name, sep=""), width=20, justify="left" )
        arg.value <- ""
        arg.value <- try(as.character( as.character(model.FUN.args[i])), TRUE)
             
        writeLines(c(arg.name, ":", arg.value), InfoTXT.TextFile, sep=" ") 
        writeLines("", InfoTXT.TextFile) # writing a blank line with a carriage return             
      } # FOR end
  } # IF end
  # Closing the text file
  close(InfoTXT.TextFile) 
  
  ########################################################################  
  #                        Text Files initialization                     #
  ########################################################################  

  # File 'Verification-ModelOut.txt' #
  model.out.text.fname <- paste(file.path(drty.out), "/", "Verification-ModelOut.txt", sep="")
  OFout.Text.file  <- file(model.out.text.fname, "w+")
  #c(isOpen(Tfile, "r"), isOpen(Tfile, "w")) # both TRUE
  #writeLines( paste("t", 1:length(obs)), OFout.Text.file, sep=" ")   
  writeLines(c("Param", "GoF", "Model_Output"), OFout.Text.file, sep="  ") 
  writeLines("", OFout.Text.file) # writing a blank line with a carriage return
  close(OFout.Text.file)  
        

  # File 'Verification-ParamValues.txt' #
  # with the parameters values for each partcile in each iteration
  gof.text.fname <- paste(file.path(drty.out), "/", "Verification-ParamValues.txt", sep="")
  Params.Text.file  <- file(gof.text.fname, "w+")           
  #c(isOpen(Tfile, "r"), isOpen(Tfile, "w")) # both TRUE
  writeLines(c("ParamNmbr", gof.name, Parameter.names), Params.Text.file, sep=" ")
  writeLines("", Params.Text.file) # writing a blank line with a carriage return
  close(Params.Text.file)          
        

  ##############################################################################
  #                            MAIN BODY
  ##############################################################################
  
  gof.all <- numeric(nparamsets)
  
  for ( p in 1:nparamsets) {
  
    if (verbose) message("                    |                      ")
    if (verbose) message("                    |                      ") 
    if (verbose) message("==============================================================")
    if (verbose) message( paste("[ Running parameter set ", 
                                format( p, width=4, justify="left" ), 
                                "/", nparamsets, " => ", 
                                format( round(100*p/nparamsets,2), width=7, justify="left" ), "%",
                                ". Starting... ]", sep="") )
    if (verbose) message("==============================================================")
  
    # Opening the file 'Verification-ModelOut.txt' for appending
    OFout.Text.file <- file(model.out.text.fname, "a")   

    # Opening the file 'Verification-ParamValues.txt' for appending
    Params.Text.file <- file(gof.text.fname, "a") 

    # Getting the parameter values  
    param.values <- as.numeric(par[p,])
    
    ############################################################################
    # 2)                       Running the hydrological model                  #
    ############################################################################
    
    # If the user wants to create a plot with obs vs sim
    if (do.plots) {
      do.png         <- TRUE
      png.fname      <- paste(drty.out, "/ParameterSet_", p, ".png", sep="")
      main <- paste("Parameter Set:", p, sep="")
      model.FUN.args <- modifyList(model.FUN.args, list(do.png=do.png, png.fname=png.fname, main=main)) 
    } # IF end
    
    # Model evaluation 
    if ( fn.name == "hydromod" ) {
      model.FUN.args <- modifyList(model.FUN.args, list(param.values=param.values)) 
      hydromod.out   <- do.call(model.FUN, as.list(model.FUN.args)) 
    } else hydromod.out <- do.call(fn, list(param.values))
     
        
    ############################################################################
    # 3)                     Extracting simulated values                       #                                 
    ############################################################################
                  
    # Extracting the simulated values and the goodness of fit
    if ( fn.name == "hydromod" ) {
      sims <- as.numeric(hydromod.out[["sim"]])
      GoF  <- as.numeric(hydromod.out[["GoF"]])
    } else {
        sims <- as.numeric(hydromod.out)
        GoF  <- as.numeric(hydromod.out)
      } # ELSe end
        
   
    gof.all[p] <- GoF

    # Writing to the 'Verification-ModelOut.txt' file
    suppressWarnings( tmp <- formatC(GoF, format="E", digits=digits, flag=" ") )
    if ( fn.name != "hydromod" ) {
      writeLines(as.character(c(p, tmp, tmp)), OFout.Text.file, sep="  ") 
    } else suppressWarnings( writeLines(as.character(c(p, tmp, formatC(sims, format="E", digits=digits, flag=" ") )), OFout.Text.file, sep="  ") )
    writeLines("", OFout.Text.file) # writing a blank line with a carriage return  
    
    # Writing to the 'Verification-ParamValues.txt' file
    suppressWarnings( writeLines(as.character(c(p, tmp, formatC(param.values, format="E", digits=digits, flag=" ") )), Params.Text.file, sep="  ") )
    writeLines("", Params.Text.file) # writing a blank line with a carriage return
    
    # Closing the output text files
    close(OFout.Text.file)
    close(Params.Text.file) 
    
  } # FOR end
 
  if (verbose) message("                              |                               ")
  if (verbose) message("                              |                               ") 
  if (verbose) message("==============================================================")
  if (verbose) message("[================ Verification finished ! ===================]")
  if (verbose) message("==============================================================")
  
  ##############################################################################
  # 4)                    Writing Ending and Elapsed Time                      #                                 
  ##############################################################################
  InfoTXT.TextFile <- file(InfoTXT.fname, "a")    
  #c(isOpen(Tfile, "r"), isOpen(Tfile, "w")) # both TRUE
  writeLines("================================================================================", InfoTXT.TextFile) # writing a separation line with a carriage return
  writeLines(c("Ending Time            :", date()), InfoTXT.TextFile, sep=" ")
  writeLines("", InfoTXT.TextFile) # writing a blank line with a carriage return
  writeLines("================================================================================", InfoTXT.TextFile) # writing a separation line with a carriage return
  # Total time of the simulations
  Time.Fin <- Sys.time()
  writeLines(c("Elapsed Time           :", format(round(Time.Fin - Time.Ini, 2))), InfoTXT.TextFile, sep=" ")
  writeLines("", InfoTXT.TextFile) # writing a blank line with a carriage return
  writeLines("================================================================================", InfoTXT.TextFile) # writing a separation line with a carriage return
  close(InfoTXT.TextFile)
  
  ##############################################################################
  # 5)                    Creating the output                                  #                                 
  ##############################################################################
  ifelse(MinMax=="min", best.rowindex <- which.min(gof.all),
                        best.rowindex <- which.max(gof.all)) 

  best.gof <- gof.all[best.rowindex]                       
  best.par <- par[best.rowindex,]

  out <- list(gofs=gof.all, best.gof=best.gof, best.par=best.par )

  return(out)

} # 'verification' END
