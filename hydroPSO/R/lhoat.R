# File lhoat.R
# Part of the hydroPSO R package, http://www.rforge.net/hydroPSO/ ; 
#                                 http://cran.r-project.org/web/packages/hydroPSO
# Copyright 2011-2014 Mauricio Zambrano-Bigiarini & Rodrigo Rojas
# Distributed under GPL 2 or later


################################################################################
###  Function for evaluating the hydrological model for a single particle  #####
################################################################################
### Started: 13-May-2013                                                     ###
### Updates: 14-May-2013 ; 16-May-2013 ; 04-Jun-2013                         ###
###          07-Feb-2014                                                     ###
################################################################################
hydromod.eval.SA <- function(j, Thetas, nparamsets,
                             N, X.Boundaries,
                             write2disk=FALSE,
                             REPORT, verbose, digits, 
                             model.FUN, model.FUN.args, 
                             parallel, ncores, part.dirs) {

  if (parallel!="none")         
     model.FUN.args <- modifyList(model.FUN.args, list(model.drty=part.dirs[j]) ) 
    
  # Creating the R output
  nelements <- 2        
  out       <- vector("list", nelements)
    
  gof.is.numeric <- FALSE
    
  while (!gof.is.numeric) {
  
    # j-th parameter set
    params       <- Thetas[j,]
    suppressWarnings( param.values <- as.numeric(formatC(params, format="E", digits=digits)) )

    ############################################################################
    # 4)                       Running the hydrological model                  #
    ############################################################################

    # Evaluating the hydrological model
    model.FUN.args <- modifyList(model.FUN.args, list(param.values=params) ) 
    hydromod.out   <- do.call(model.FUN, as.list(model.FUN.args)) 
    
#    # If the user wants to create a plot with obs vs sim
#    if (do.plots) {
#      do.png         <- TRUE
#      png.fname      <- paste(drty.out, "/LHS-Point_0_", j, ".png", sep="")
#      k              <- ceiling(P/5)
#      title          <- character(k)
#      for (m in 1:k) {
#        namess   <- format(names(params[(5*(m-1)+1):(5*m)]), 10, justify="left")
#        values   <- format(round(params[(5*(m-1)+1):(5*m)], 3), 7, justify="left")
#        title[m] <- paste(namess, "=", values, " ; ", collapse="")
#      } # FOR end
#      main <- paste(title, sep="\n")
#      model.FUN.args <- modifyList(model.FUN.args, list(do.png=do.png, png.fname=png.fname, main=main)) 
#    } # IF end

        
    ############################################################################
    # 5)                     Extracting simulated values                       #                                 
    ############################################################################
                  
    out[[1]] <- as.numeric(hydromod.out[["GoF"]])
    out[[2]] <- hydromod.out[["sim"]]
        
    # Finding out if the GoF had a finite value or not
    ifelse(is.finite(out[[1]]), gof.is.numeric <- TRUE, gof.is.numeric <- FALSE)
    
    # If the current point of the initial LHS leads to a GoF=NA, it is replaced
    if (!gof.is.numeric) {
      warning(" parameter set ", j, ": not numeric GoF ! => it was replaced")
      tmp        <- rLHS(n=N, ranges=X.Boundaries)
      Thetas[j,] <- tmp[j,]
    } # IF end
        
  } # WHILE '(!gof.is.numeric)' end
  
  # meaningful names
  names(out)[1:nelements] <- c("GoF", "model.out") 

  if ( j/REPORT == floor(j/REPORT) ) {
    if (verbose) message( "[ Parameter set ", 
                                format( j, width=4, justify="left" ), 
                                "/", nparamsets, 
                          ". Finished.   GoF: ", format(hydromod.out[["GoF"]], scientific=TRUE, digits=digits), 
                          "]" )  
  } # IF end
        
  return(out)
  
} # 'hydromod.eval.SA' END


################################################################################
#                  Latin-Hypercube One-At-a-Time                               #
################################################################################
# Purpose   : Run a sensitivity analysis for the parameters of the             #
#             hydrological model, by using Latin-Hypercupe One-factor-At-a-Time#  
#             developed by  van Griensven et al., 2006                         #
################################################################################
# Reference : A. van Griensven, T. Meixner, S. Grunwald, T. Bishop, M. Diluzio, 
#             R. Srinivasan, A global sensitivity analysis tool for the 
#             parameters of multi-variable catchment models, Journal of Hydrology, 
#             Volume 324, Issues 1-4, 15 June 2006, Pages 10-23, ISSN 0022-1694, 
#             DOI: 10.1016/j.jhydrol.2005.09.008.
#             (http://www.sciencedirect.com/science/article/pii/S0022169405004488)
################################################################################
# Output: A list of two elements:                                             #
#        1) ParameterSets: A matrix with all the parameter sets used in the    #
#                          LH-OAT                                              #
#        2) Ranking      : data.frame with 4 columns:                          #
#           2.1) RankingNmbr       : integer, with the ranking of the sensitive#
#                                    parameters, 1 for the most sensitive para-#
#                                    meter                                     #  
#           2.2) ParameterName     : character with the name of each one of the#
#                                    parameters, sorted in decreasing order,   #
#                                    from the most sensitive to the least one. #
#           2.3) RelativeImportance: numeric with the relative importance of   #
#                                    each one of the parameters, sorted in     #
#                                    decreasing order, from the most sensitive #
#                                    to the least one.                         #
#           2.4) RelativeImportance.Norm: numeric with the normalised relative #
#                                    importance of each one of the parameters, #
#                                    sorted in decreasing order, from the most #
#                                    sensitive to the least one.               #
################################################################################
# Author  : Mauricio Zambrano-Bigiarini                                        #
# Started : 23-Jun-2011                                                        #
# Updates : 26-Jan-2012 ; 02-Feb-2012 ; 13-Feb-2012 ; 23-Feb-2012              #
#           09-May-2013 ; 13-May-2013 ; 15-May-2013 ; 16-May-2013              #
#           28-May-2013 ; 27-Aug-2013 ; 27-Dec-2013                            #
#           07-Feb-2014 ; 09-Abr-2014                                          #
################################################################################

lhoat <- function(
                  fn="hydromod",  
                  ...,  
                  lower=-Inf,
                  upper=Inf,
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
  } else 
      if ( is.character(fn) | is.function(fn) )  {
        if (is.character(fn)) {
          fn.name <- fn
          fn      <- match.fun(fn)
	} else if (is.function(fn)) {
	    fn.name <- as.character(substitute(fn))
	    fn      <- fn
	  } # ELSE end
      } else stop("Missing argument: 'class(fn)' must be in c('function', 'character')")      
        
  # checking length of 'lower' and 'upper'
  if (length(lower) != length(upper) )
    stop( paste( "Invalid argument: 'length(lower) != length(upper) (", length(lower), "!=", length(upper), ")'", sep="" ) )
        
          
  ########################################################################        
  con <- list(
        
          N=3,                            # number of strata to be used for each parameter
          f=0.15,                         # fraction by which each single parameter is changed within the Morris OAT design
                
          drty.in="PSO.in",
          drty.out="LH_OAT",              # Character, with the name of the directory that will store the results of the LH-OAT. 
          param.ranges="ParamRanges.txt", # Character, with the name of the file that stores the desired range of variation for each parameter                          
          digits=7,
          normalise=FALSE,     
          normaliseRanking=FALSE,

          gof.name="GoF",
          do.plots=FALSE,
          write2disk=TRUE,
          verbose= TRUE,                   # logical, indicating if progress messages have to be printed          
          REPORT=100, 
	  parallel=c("none", "multicore", "parallel", "parallelWin"),
	  par.nnodes=NA,
	  par.pkgs= c()
             )
              
  parallel <- match.arg(control[["parallel"]], con[["parallel"]]) 
  
  nmsC <- names(con)
 
  con[(namc <- names(control))] <- control
  
  if (length(noNms <- namc[!namc %in% nmsC])) 
    warning("[Unknown names in control: ", paste(noNms, collapse = ", "), " (not used) !]")
    
  N                <- con[["N"]]
  f                <- con[["f"]]
  drty.in          <- con[["drty.in"]]
  drty.out         <- con[["drty.out"]]
  param.ranges     <- con[["param.ranges"]]         
  digits           <- con[["digits"]]
  normalise        <- as.logical(con[["normalise"]])   

  gof.name       <- con[["gof.name"]]
  do.plots       <- as.logical(con[["do.plots"]])  
  write2disk     <- as.logical(con[["write2disk"]])
  verbose        <- as.logical(con[["verbose"]])
  REPORT         <- con[["REPORT"]] 
  par.nnodes     <- con[["par.nnodes"]]
  par.pkgs       <- con[["par.pkgs"]] 
        
  ########################################################################
  ##################### Dummy checkings ##################################

  # Checking that 'N' is integer
  if ( trunc(N) != N ) stop( "Invalid argument: 'N' must be integer !" )
  
  # Checking that '0 < f < 1' 
  if ( (f <= 0) | (f >= 1) ) stop( "Invalid argument: 'f' must be in ]0, 1[" )   
        
  # 'hydromod' checkings
  if (fn.name=="hydromod") {

    # Checking that 'param.ranges' really exists
    if ( !file.exists( param.ranges ) )
       stop( paste("Invalid argument: The file '", param.ranges, "' does not exist !", sep="") ) 
             
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


  # checking 'X.Boundaries' 
  if (fn.name=="hydromod") {
        
     if (verbose) message("==============================================================")
     if (verbose) message("[   1)   Reading 'param.ranges' ...                          ]")
     if (verbose) message("==============================================================")
  
      # Reading the desired range of variation for each parameter 
      X.Boundaries <- read.ParameterRanges(ParamRanges.fname=param.ranges) 

      lower <- X.Boundaries[,1]
      upper <- X.Boundaries[,2]

          
  } else {
      if ( (lower[1L] == -Inf) || (upper[1L] == Inf) ) {
      #if (any(upper==Inf | lower==-Inf))
        stop( "Invalid argument: 'lower' and 'upper' boundaries must be finite !!'" )
      } else X.Boundaries <- cbind(lower, upper)              
    } # ELSE end
          
  # Computing 'P', the Dimension of the Solution Space
  P <- nrow(X.Boundaries)

  # Meaningful name of each one of the parameters
  if (is.null(rownames(X.Boundaries))) {
    param.IDs <- paste("Param", 1:P, sep="")
  } else param.IDs <- rownames(X.Boundaries)   
  
  # Total Number of parameter sets to be used in the LH-OAT
  nparamsets  <- (P+1)*N

  # Checking report
  if (nparamsets < REPORT) {
      REPORT <- nparamsets
      warning("[ 'REPORT' is greater than 'nparamsets' => 'REPORT=nparamsets' ]")
  } # IF end 
  
  if (verbose) message("                                                              ")
  if (verbose) message("[ Number of strata for each parameter (N) : ", N, " ]")    
  
  if (verbose) message("                                                              ")
  if (verbose) message("[ Number of Parameter Sets to be run      : ", nparamsets, " ]") 

  if (normalise) {
      # Backing up the original boundaries
      lower.ini <- lower
      upper.ini <- upper
      X.Boundaries.ini <- X.Boundaries
      LOWER.ini <- matrix( rep(lower.ini, nparamsets), nrow=nparamsets, byrow=TRUE)
      UPPER.ini <- matrix( rep(upper.ini, nparamsets), nrow=nparamsets, byrow=TRUE)
      
      # normalising
      lower <- rep(0, P)
      upper <- rep(1, P)
      X.Boundaries <- cbind(lower, upper)
      rownames(X.Boundaries) <- param.IDs
    } # IF end
        
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
  ##                                  parallel                                 #
  ##############################################################################
  if (parallel != "none") {
    
    if ( ( (parallel=="multicore") | (parallel=="parallel") ) & 
       ( (R.version$os=="mingw32") | (R.version$os=="mingw64") ) )
       stop("[ Fork clusters are not supported on Windows =>  'parallel' can not be set to '", parallel, "' ]")
    
    ifelse(parallel=="parallelWin", parallel.pkg <- "parallel",  parallel.pkg <- parallel)                
    if ( !require(parallel) ) {
            warning("[ Package '", parallel.pkg, "' is not installed =>  parallel='none' ]")
            parallel <- "none"
    }  else { 
      
         if (verbose) message("                               ")
         if (verbose) message("[ Parallel initialization ... ]")
      
         fn1 <- function(i, x) fn(x[i,])
      
         nnodes.pc <- detectCores()
         if (verbose) message("[ Number of cores/nodes detected: ", nnodes.pc, " ]")
      
         if ( (parallel=="parallel") | (parallel=="parallelWin") ) {               
            logfile.fname <- paste(file.path(drty.out), "/", "parallel_logfile.txt", sep="") 
            if (file.exists(logfile.fname)) file.remove(logfile.fname)
         } # IF end
                      
         if (is.na(par.nnodes)) {
           par.nnodes <- nnodes.pc
         } else if (par.nnodes > nnodes.pc) {
               warning("[ 'nnodes' > number of detected cores (", par.nnodes, ">", nnodes.pc, ") =>  par.nnodes=", nnodes.pc, " ] !",)
               par.nnodes <- nnodes.pc
           } # ELSE end
         if (par.nnodes > nparamsets) {
           warning("[ 'par.nnodes' > N*(P+1) (", par.nnodes, ">", nparamsets, ") =>  par.nnodes=", nparamsets, " ] !")
           par.nnodes <- nparamsets
         } # ELSE end  
         if (verbose) message("[ Number of cores/nodes used    : ", par.nnodes, " ]")                 
              
         if (parallel=="parallel") {
             ifelse(write2disk, 
                    cl <- makeForkCluster(nnodes = par.nnodes, outfile=logfile.fname),
                    cl <- makeForkCluster(nnodes = par.nnodes) )         
         } else if (parallel=="parallelWin") {      
             ifelse(write2disk,
                 cl <- makeCluster(par.nnodes, outfile=logfile.fname),
                 cl <- makeCluster(par.nnodes) )
             pckgFn <- function(packages) {
               for(i in packages) library(i, character.only = TRUE)
             } # 'packFn' END
             clusterCall(cl, pckgFn, par.pkgs)
             clusterExport(cl, ls.str(mode="function",envir=.GlobalEnv) )
             if (fn.name=="hydromod") {
               clusterExport(cl, model.FUN.args$out.FUN)
               clusterExport(cl, model.FUN.args$gof.FUN)
             } # IF end                   
           } # ELSE end                   
                            
         if (fn.name=="hydromod") {
           if (!("model.drty" %in% names(formals(hydromod)) )) {
              stop("[ Invalid argument: 'model.drty' has to be an argument of the 'hydromod' function! ]")
           } else { # Copying the files in 'model.drty' as many times as the number of cores
           
               model.drty <- path.expand(model.FUN.args$model.drty)
                 
               files <- list.files(model.drty, full.names=TRUE, include.dirs=TRUE) 
               tmp <- which(basename(files)=="parallel")
               if (length(tmp) > 0) files <- files[-tmp]
               parallel.drty <- paste(file.path(model.drty), "/parallel", sep="")

               if (file.exists(parallel.drty)) {                      
                 if (verbose) message("[ Removing the 'parallel' directory ... ]")    
                 try(unlink(parallel.drty, recursive=TRUE, force=TRUE))
               } # IF end 
               dir.create(parallel.drty)

               mc.dirs <- character(par.nnodes)
               if (verbose) message("                                                     ")
               for (i in 1:par.nnodes) {
                 mc.dirs[i] <- paste(parallel.drty, "/", i, "/", sep="")
                 dir.create(mc.dirs[i])
                 if (verbose) message("[ Copying model input files to directory '", mc.dirs[i], "' ... ]")
                 file.copy(from=files, to=mc.dirs[i], overwrite=TRUE, recursive=TRUE)
               } # FOR end
                 
               n         <- ceiling(nparamsets/par.nnodes)        
               part.dirs <- rep(mc.dirs, n)[1:nparamsets]  
             } # ELSE end                 
         } # IF end
           
       } # ELSE end  
  
  }  # IF end    
  ##############################################################################
    
  
  ##############################################################################
  #                            Writing Info File
  ##############################################################################  
  if (write2disk) {   
    # File 'Verification-logfile.txt' #        
    InfoTXT.fname <- paste(file.path(drty.out), "/", "LH_OAT-logfile.txt", sep="")
    InfoTXT.TextFile  <- file(InfoTXT.fname , "w+")
    #c(isOpen(Tfile, "r"), isOpen(Tfile, "w")) # both TRUE
    writeLines("================================================================================", InfoTXT.TextFile) # writing a separation line with a carriage return
    writeLines(c("hydroPSO version     :", sessionInfo()$otherPkgs$hydroPSO$Version), InfoTXT.TextFile, sep="  ")
    writeLines("", InfoTXT.TextFile) #writing a blank line with a carriage return
    writeLines(c("hydroPSO Built       :", sessionInfo()$otherPkgs$hydroPSO$Built), InfoTXT.TextFile, sep="  ")
    writeLines("", InfoTXT.TextFile) # writing a blank line with a carriage return
    writeLines(c("R version            :", sessionInfo()[[1]]$version.string), InfoTXT.TextFile, sep="  ")
    writeLines("", InfoTXT.TextFile) # writing a blank line with a carriage return
    writeLines(c("Platform             :", sessionInfo()[[1]]$platform), InfoTXT.TextFile, sep="  ")
    writeLines("", InfoTXT.TextFile) # writing a blank line with a carriage return
    writeLines("================================================================================", InfoTXT.TextFile) # writing a separation line with a carriage return
    Time.Ini <- Sys.time()
    writeLines(c("Starting Time        :", date()), InfoTXT.TextFile, sep="  ")
    writeLines("", InfoTXT.TextFile) # writing a blank line with a carriage return
    writeLines("================================================================================", InfoTXT.TextFile) # writing a separation line with a carriage return
    writeLines(c("N (number of strata) :", N), InfoTXT.TextFile, sep=" ") 
    writeLines("", InfoTXT.TextFile) # writing a blank line with a carriage return
    writeLines(c("f (changing factor)  :", f), InfoTXT.TextFile, sep=" ") 
    writeLines("", InfoTXT.TextFile) # writing a blank line with a carriage return
    writeLines(c("drty.in              :", drty.in), InfoTXT.TextFile, sep=" ") 
    writeLines("", InfoTXT.TextFile) # writing a blank line with a carriage return
    writeLines(c("drty.out             :", drty.out), InfoTXT.TextFile, sep=" ") 
    writeLines("", InfoTXT.TextFile) # writing a blank line with a carriage return
    writeLines(c("param.ranges         :", param.ranges), InfoTXT.TextFile, sep=" ") 
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
    writeLines(c("normalise            :", normalise), InfoTXT.TextFile, sep=" ") 
    writeLines("", InfoTXT.TextFile) # writing a blank line with a carriage return
    writeLines(c("parallel             :", parallel), InfoTXT.TextFile, sep=" ")  
    writeLines("", InfoTXT.TextFile)  
    if (parallel!="none") {
      writeLines(c("par.nnodes           :", par.nnodes), InfoTXT.TextFile, sep=" ") 
      writeLines("", InfoTXT.TextFile)
      writeLines(c("par.pkgs             :", par.pkgs), InfoTXT.TextFile, sep=" ") 
      writeLines("", InfoTXT.TextFile)     
    } # IF end
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
  } # IF end       

  
  ########################################################################  
  #                        Text Files initialization                     #
  ########################################################################  
  if (write2disk) { 
    # File 'LH_OAT-out.txt' #
    model.out.text.fname <- paste(file.path(drty.out), "/", "LH_OAT-out.txt", sep="")
    model.out.text.file  <- file(model.out.text.fname, "w+")
    #c(isOpen(Tfile, "r"), isOpen(Tfile, "w")) # both TRUE
    #writeLines( paste("t", 1:length(obs)), model.out.text.file, sep=" ") 
    close(model.out.text.file)  
        

    # File 'LH_OAT-gof.txt' #
    # with the parameters values for each partcile in each iteration
    gof.text.fname <- paste(file.path(drty.out), "/", "LH_OAT-gof.txt", sep="")
    gof.text.file  <- file(gof.text.fname, "w+")           
    #c(isOpen(Tfile, "r"), isOpen(Tfile, "w")) # both TRUE
    writeLines(c(gof.name, param.IDs), gof.text.file, sep=" ") 
    writeLines("", gof.text.file) # writing a blank line with a carriage return
    close(gof.text.file)          
  } # IF end
        

  ##############################################################################
  #                            MAIN BODY
  ##############################################################################
  
  if (verbose) message("                                                              ")
  if (verbose) message("==============================================================")
  if (verbose) message("[  2)   Initial  LHS ...                                     ]")
  if (verbose) message("==============================================================")
  
  # Initial N points with LHS
  Theta.Ini <- rLHS(n=N, ranges=X.Boundaries)
  
  if (verbose) message("                                                              ")
  if (verbose) message("==============================================================")
  if (verbose) message("[  3)   Running LH-OAT ...                                   ]")
  if (verbose) message("==============================================================")
  
  # Output Sensitivity Matrix
  S           <- matrix(NA, ncol=P, nrow=N)
  colnames(S) <- param.IDs
  
  # Parameter sets that will be used in the LH-OAT
  Thetas           <- matrix(NA, nrow=nparamsets, ncol=P)
  colnames(Thetas) <- param.IDs

  Theta.index <- seq(1, nparamsets, by=(P+1))
  for (j in 1:N) {
    Thetas[Theta.index[j], ] <- Theta.Ini[j, ]
    for (i in 1:P) {
      canonical    <- rep(1, P)
      canonical[i] <- 1 + f*sign(rnorm(1))
      Theta.New    <- Theta.Ini[j,]*canonical
      #n            <- n+1
      Thetas[Theta.index[j] + i, ] <- Theta.New   
      #param.values <- as.numeric(formatC(Theta.New, format="E", digits=digits))
    } # FOR 'i' end
  } # FOR 'j' end

  # Goodness-of-fit of each paremter set used in the LH-OAT
  gof    <- rep(NA, nparamsets)

  ModelOut <- vector("list", nparamsets)
  
  ##############################################################################
  #                4) Loop for each parameter set                              #
  ##############################################################################
  if (write2disk) { 
    # Opening the file 'LH_OAT-out.txt' for appending
    model.out.text.file <- file(model.out.text.fname, "a")   
    # Opening the file 'LH_OAT-gof.txt' for appending
    gof.text.file <- file(gof.text.fname, "a") 
  } else {
          model.out.text.file <- ""
          gof.text.file       <- ""
         } # ELSE end
  
  if (normalise) {
    Xn <- Thetas * (UPPER.ini - LOWER.ini) + LOWER.ini
  } else Xn <- Thetas
      
  # 3.a) Evaluate the particles fitness
  if ( fn.name != "hydromod" ) {
         
     # Evaluating an R Function 
     if (parallel=="none") {
       GoF <- apply(Xn, fn, MARGIN=1, ...)
     } else             
        if (parallel=="multicore") {
          GoF <- unlist(mclapply(1:nparamsets, FUN=fn1, x=Xn, ..., mc.cores=par.nnodes, mc.silent=TRUE)) 
        } else if ( (parallel=="parallel") | (parallel=="parallelWin") ) {
            GoF <- parRapply(cl= cl, x=Xn, FUN=fn, ...)
          } # ELSE end
	 
     gof[1:nparamsets]      <- GoF
     ModelOut[1:nparamsets] <- GoF  ###

  } else { # fn.name = "hydromod"       
	     
     if (parallel=="none") {
           out <- lapply(1:nparamsets, hydromod.eval.SA,      
                         Thetas=Xn, 
                         nparamsets=nparamsets, 
                         N=N, X.Boundaries=X.Boundaries,
                         write2disk=write2disk,
                         REPORT=REPORT, 
                         verbose=verbose, 
                         digits=digits, 
                         model.FUN=model.FUN, 
                         model.FUN.args=model.FUN.args, 
                         parallel=parallel, 
                         ncores=par.nnodes, 
                         part.dirs=mc.dirs  
                         )
                   
     } else if ( (parallel=="parallel") | (parallel=="parallelWin") ) {
                 
              out <- clusterApply(cl=cl, x=1:nparamsets, fun= hydromod.eval.SA,                                  
                                        Thetas=Xn, 
                                        nparamsets=nparamsets, 
                                        N=N, X.Boundaries=X.Boundaries,
                                        write2disk=write2disk,
                                        REPORT=REPORT, 
                                        verbose=verbose, 
                                        digits=digits, 
                                        model.FUN=model.FUN, 
                                        model.FUN.args=model.FUN.args, 
                                        parallel=parallel, 
                                        ncores=par.nnodes, 
                                        part.dirs=part.dirs                          
                                        ) # sapply END                 
                                  
             } else if (parallel=="multicore") {
                   
                       out <- mclapply(1:nparamsets, hydromod.eval.SA,       
                                                  Thetas=Xn, 
                                                  nparamsets=nparamsets, 
                                                  N=N, X.Boundaries=X.Boundaries,
                                                  write2disk=write2disk,
                                                  REPORT=REPORT, 
                                                  verbose=verbose, 
                                                  digits=digits, 
                                                  model.FUN=model.FUN, 
                                                  model.FUN.args=model.FUN.args, 
                                                  parallel=parallel, 
                                                  ncores=par.nnodes, 
                                                  part.dirs=part.dirs,
                                                  mc.cores=par.nnodes,
                                                  mc.silent=TRUE,
                                                  mc.cleanup=TRUE                     
                                                  ) # mclapply END
                                      
                     } # ELSE end
                                       
             for (j in 1:nparamsets){         
                   gof[j]              <- out[[j]][["GoF"]] 
                   ModelOut[[j]]       <- out[[j]][["model.out"]]  
                   #nfn <- nfn + 1 
                   #if(is.finite(GoF)) nfn.eff <- nfn.eff + 1  
  
                   if (write2disk) { 
                     # j-th parameter set
                     suppressWarnings( param.values <- as.numeric(formatC(Xn[j,], format="E", digits=digits)) )

                     # Writing to the 'LH_OAT-out.txt' file
                     writeLines(as.character(out[[j]][["model.out"]]), model.out.text.file, sep=" ") 
                     writeLines("", model.out.text.file) # writing a blank line with a carriage return
                     flush(model.out.text.file) 
    
                     # Writing to the 'LH_OAT-gof.txt' file
                     suppressWarnings(
                     writeLines( as.character( c(formatC(out[[j]][["GoF"]], format="E", digits=digits, flag=" "), # GoF
                                                 formatC(param.values, format="E", digits=digits, flag=" ")                                            
                                                ) ), gof.text.file, sep="  ") 
                     )
                                             
                     writeLines("", gof.text.file) # writing a blank line with a carriage return
                     flush(gof.text.file) 
                   } # IF end    
                   
             } #FOR part end               

	} # ELSE end


  ##############################################################################
  # 5)                    Writing output files                                 #
  ##############################################################################
  if (write2disk) { 
    
    if ( (fn.name != "hydromod") ) {
      if (verbose) message("==============================================================")
      if (verbose) message("[  5)   Writing output files ...                             ]")
      if (verbose) message("==============================================================")

      for (j in 1:nparamsets) {
        # Writing to the 'LH_OAT-out.txt' file
        writeLines(as.character(gof[j]), model.out.text.file, sep=" ") 
        writeLines("", model.out.text.file) # writing a blank line with a carriage return
        flush(model.out.text.file) 

        # Writing to the 'LH_OAT-gof.txt' file
        if (normalise) {
          temp <- Thetas[j,] * (upper.ini - lower.ini) + lower.ini
        } else temp <- Thetas[j,]
        temp.gof <- gof[j]
	if(is.finite(temp.gof)) {
          suppressWarnings(
	  writeLines( as.character( c(formatC(temp.gof, format="E", digits=digits, flag=" "), 
	                              formatC(temp, format="E", digits=digits, flag=" ")	                                                            
	                          ) ), gof.text.file, sep="  ") 
          )
	} else suppressWarnings( writeLines( as.character( c("NA",
	                                   formatC(temp, format="E", digits=digits, flag=" ")                                                                                  
	                               ) ), gof.text.file, sep="  ")   
                               )                                          
        writeLines("", gof.text.file) # writing a blank line with a carriage return
        flush(gof.text.file) 
      } # FOR end  
    } # IF end

    # Closing the output text files
    close(model.out.text.file)
    close(gof.text.file) 
  } # IF end

  ##############################################################################
  # 6)                    Updating the sensitivity matrix                      #                                 
  ##############################################################################

  for (j in 1: N) {
    M.Zero <- gof[Theta.index[j]]
    for (i in 1:P) {
      M.New <- gof[Theta.index[j] + i]
      S[j,i] <- abs( (M.New - M.Zero) / ( (M.New + M.Zero) / 2) ) * ( 100 / f)
    } # FOR 'i' end
  } # FOR 'j' end
  
  ##############################################################################
  # 7)                    Sensitivity of each Parameter                        #                                 
  ##############################################################################
  
  # Mean value for each parameter
  Ranking <- colMeans(S, na.rm=TRUE)
  
  if (verbose) message("                             |                                ")
  if (verbose) message("                             |                                ") 
  if (verbose) message("==============================================================")
  if (verbose) message("[==================    LH-OAT finished !    =================]")
  if (verbose) message("==============================================================")
  
  
  ##############################################################################
  ##                                   parallel                                #
  ##############################################################################
  if (parallel!="none") {
    if ( (parallel=="parallel") | (parallel=="parallelWin") )   
         stopCluster(cl)   
    if (fn.name=="hydromod") {
      if (verbose) message("                                         ")
      if (verbose) message("[ Removing the 'parallel' directory ... ]")    
      unlink(dirname(mc.dirs[1]), recursive=TRUE)
    } # IF end
  } # IF end
  
  
  ##############################################################################
  # 8)                    Writing Ending and Elapsed Time                     #                                 
  ##############################################################################
  if (write2disk) { 
    InfoTXT.TextFile <- file(InfoTXT.fname, "a")    
    #c(isOpen(Tfile, "r"), isOpen(Tfile, "w")) # both TRUE
    writeLines("================================================================================", InfoTXT.TextFile) # writing a separation line with a carriage return
    writeLines(c("Ending Time          :", date()), InfoTXT.TextFile, sep=" ")
    writeLines("", InfoTXT.TextFile) # writing a blank line with a carriage return
    writeLines("================================================================================", InfoTXT.TextFile) # writing a separation line with a carriage return
    # Total time of the simulations
    Time.Fin <- Sys.time()
    writeLines(c("Elapsed Time         :", format(round(Time.Fin - Time.Ini, 2))), InfoTXT.TextFile, sep=" ")
    writeLines("", InfoTXT.TextFile) # writing a blank line with a carriage return
    writeLines("================================================================================", InfoTXT.TextFile) # writing a separation line with a carriage return
    close(InfoTXT.TextFile)
  } # IF end
  
  ##############################################################################
  # 9)        Computing the Final Ranking of Sensitivity                       #                                 
  ##############################################################################
  # Sorting the parameters, from the most sensitive to the least one
  Ranking <- sort(Ranking, decreasing=TRUE, na.last=TRUE)
  
  # Adding a column with a sequential number for the ranking
  Ranking <- data.frame(RankingNmbr=format(as.character(1:P), width=11, justify="left"), 
                        ParameterName=format(names(Ranking), width=13, justify="left"), 
                        RelativeImportance=as.numeric(Ranking) )                        
  Ranking[, "RankingNmbr"] <- as.character(Ranking[, "RankingNmbr"])
  Ranking[,"RelativeImportance.Norm"] <- Ranking[,"RelativeImportance"]/sum(Ranking[,"RelativeImportance"], na.rm=TRUE)
  
                        
  # Assigning the same worst ranking to all the insensitive parameters
  row.index <- which(Ranking[,"RelativeImportance"]==0)
  ninsens   <- length(row.index)
  if (ninsens > 0)
    Ranking[row.index, "RankingNmbr"] <- format(as.character(rep(P, ninsens)), width=11, justify="left")
  
  ##############################################################################
  # 10)                    Creating the output                                 #                                 
  ##############################################################################
  
  if (write2disk) { 
    # Writing the output file 'LH_OAT-Ranking.txt'
    Ranking.fname <- paste(file.path(drty.out), "/", "LH_OAT-Ranking.txt", sep="")
    write.table(Ranking, file=Ranking.fname, row.names=FALSE, quote=FALSE)
  } # IF end
  
  ## "pre-allocate" an empty list of length 2
  out <- vector("list", 2)

  if (normalise) {
    Xn <- Thetas * (UPPER.ini - LOWER.ini) + LOWER.ini
  } else Xn <- Thetas
  
  out[[1]] <- Xn
  out[[2]] <- Ranking
  names(out) <- c("ParameterSets", "Ranking")                    
  
  return(out)

} # 'lhoat' END
