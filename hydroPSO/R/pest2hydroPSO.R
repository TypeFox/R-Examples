# File pest2hydroPSO.R
# Part of the hydroPSO package, http://www.rforge.net/hydroPSO/
# Copyright 2012-2013 Mauricio Zambrano-Bigiarini & Rodrigo Rojas
# Distributed under GPL 2 or later

################################################################################
##                             .pst2paramranges                               ##
################################################################################
# Author : Mauricio Zambrano-Bigiarini & Rodrigo Rojas                        ##
################################################################################
# Created: 08-Nov-2012                                                        ##
# Updates: 09-Nov-2012                                                        ##
################################################################################
# Purpose  : To write the 'ParamRanges.txt' hydroPSO input file               ##
################################################################################

# 'drty.model': character, with the name of the directory with the input and exe
#               files of the model
# 'names'     : character vector, with the names of the parameters to be calibrated
# 'ini'       : numeric, with the initial guess used in PEST for each parameter
#               to be calibrated (not used in hydroPSO)
# 'min'       : numeric, with the lowest boundaries for the the parameters to be 
#               calibrated
# 'max'       : numeric, with the highest boundaries for the the parameters to be 
#               calibrated
# 'fname.out' : character, with the name of the output text file with the boudaries
#               for each parameter to be calibrated


.pst2paramranges <- function(drty.model, names, ini, min, max, 
                             fname.out="ParamRanges.txt") {

  drty.bak <- getwd()
  setwd(drty.model)
  
  nparam <- length(names)
  
  if ( (nparam!=length(ini)) | (nparam!=length(min)) | (nparam!=length(max)) )
    stop("Invalid argument: dimensions do not match !")
    
  field.names <- c("ParameterNmbr", "ParameterName", "MinValue", "MaxValue", "IniPEST")
  
  tmp <- matrix(NA, ncol=5, nrow=nparam)
  tmp <- as.data.frame(tmp)
  tmp[,1] <- 1:nparam
  tmp[,2] <- names
  tmp[,3] <- min
  tmp[,4] <- max
  tmp[,5] <- ini
  
  colnames(tmp) <- field.names
  
  write.table(tmp, file=fname.out, row.names=FALSE, quote=FALSE)
  
  setwd(drty.bak)

} # '.pst2paramranges' END



################################################################################
##                           .pst2paramfiles                                  ##
################################################################################
# Author : Mauricio Zambrano-Bigiarini & Rodrigo Rojas                        ##
################################################################################
# Created: 08-Nov-2012                                                        ##
# Updates: 09-Nov-2012                                                        ##
#          04-Jun-2013 : 05-Jun-2013                                          ##
################################################################################
# Purpose  : To write the 'ParamFiles.txt' hydroPSO input file                ##
################################################################################
.pst2paramfiles <- function(drty.model, tpls, inputs, param.names, 
                            fname.out="ParamFiles.txt", decimals=5) {

  drty.bak <- getwd()
  setwd(drty.model)
  
  ntpl    <- length(tpls)
  ninputs <- length(inputs)
  nparam  <- length(param.names)
  
  if ( (ntpl!=ninputs) )
    stop("Invalid argument: dimensions do not match !")
    
  field.names <- c("ParameterNmbr", "ParameterName", "Filename", 
                   "Row.Number", "Col.Start", "Col.End", "decimals")
  
  #tmp <- matrix(NA, ncol=length(filed.names), nrow=nparam)
  
  # output creation
  tmp <- matrix(NA, ncol=length(field.names), nrow=1)
  tmp <- as.data.frame(tmp)
  colnames(tmp) <- field.names
  
  out <- vector("list", nparam)
  for (i in 1:nparam) out[[i]] <- tmp
  
  # text width of user for showing the number of parameter on the screen
  nparam.width <- nchar(nparam) 
  
  # text width of user for showing the number of .tpl on the screen
  ntpl.width <- nchar(ntpl)   
    
  for (p in 1:nparam) {
  
    message("[ Processing parameter ",format( p, width=nparam.width, justify="left" ), "/", 
              format( nparam, width=nparam.width, justify="left" ), ": ", param.names[p], " ... ]")    
  
    for (f in 1:ntpl) {
        
      message("   [ Processing   file    ",format( f, width=ntpl.width, justify="left" ), "/", 
              format( ntpl, width=ntpl.width, justify="left" ), ": ", basename(tpls[f]), " ... ]")
  
      if (!file.exists(tpls[f])) stop("Invalid argument: file '", tpls[f], " does not exist !!")
      x <- readLines(tpls[f])
    
      # getting the token
      token <- strsplit(x[1], " ", useBytes=TRUE)[[1]][2]
    
        for (l in 2:length(x)) {          
          # Does 'param.names[p]' exist in 'x[l]' ?
          exists <- grep(param.names[p], x[l], useBytes=TRUE)             
          if (length(exists) > 0) {
          
            # Does 'token' exist in 'x[l]' ?
            token.pos <- which(strsplit(x[l], '', useBytes=TRUE)[[1]]==token)            
            if (length(token.pos)>0) {
              out[[p]] <- rbind(out[[p]], c(p, param.names[p], inputs[f], l-1, token.pos[1], token.pos[2], decimals) )
              if (length(token.pos) >2) {
                for ( t in seq(3, length(token.pos), by=2) )
                out[[p]] <- rbind(out[[p]], c(p, param.names[p], inputs[f], l-1, token.pos[t], token.pos[t+1], decimals) )
              } # IF end
            } # IF end            
          } # IF end    
        } # FOR l end
    
    } # FOR f end
    
    # removing dummy 1st row
    out[[p]] <- out[[p]][-1,]
   } # FOR 'p' end
   
  for (p in 1:nparam) tmp <- rbind(tmp, out[[p]])
  
  # removing dummy 1st row
  tmp <- tmp[-1, ]
  
  # identifying possible errors
  error.index <- which ( (as.numeric(tmp[,6]) - as.numeric(tmp[,5]) + 1 ) <=  decimals)
  if (length(error.index) > 0) {
    warning("In ParamFiles.txt, decimal places have to be manually corrected:", paste(error.index) )
  } # IF

  # Writing the output file
  write.table(tmp, file=fname.out, row.names=FALSE, quote=FALSE)
  
  setwd(drty.bak)

} # '.pst2paramfiles' END


################################################################################
##                             pest2hydroPSO                                  ##
################################################################################
# Author : Mauricio Zambrano-Bigiarini & Rodrigo Rojas                        ##
################################################################################
# Created: 08-Nov-2012                                                        ##
# Updates: 09-Nov-2012 ; 15-Nov-2012                                          ##
#          28-May-2013 ; 04-Jun-2013                                          ##
################################################################################
# Purpose  : To import the PEST input files (.pst, .tpl) to be used within    ##
#            hydroPSO (ParamFiles.txt, ParamRanges.txt, hydroPSO_Rscript.R)   ##
################################################################################
# 'pst.fname'         : character, with name of the PEST input file (.pst),   ##
#                       which contains all the information regarding          ##
#                       parameters, observations and template files (.tpl and ##
#                       .ins) used by PEST                                    ##
# 'drty.pest'         : character, path to the executable file of PEST. ALL   ##
#                       the input files required to run PEST with the model   ##
#                       have to be located within this directory (.tpl and    ##
#                       .ins)                                                 ##
# 'drty.model'        : character, path to the executable file of the model   ##
#                       specified in \code{exe.fname}. ALL the input files    ##
#                       required to run the model have to be located within   ##
#                       this directory                                        ##
# 'drty.out'          : character, name of the directory that will store all  ##
#                       the output files produced by this function            ##
# 'param.files'       : character, file name (full path) of the hydroPSO      ##
#                       input file storing the location and names of the model##
#                       files that have to be modified for each parameter     ##
#                       subject to calibration.                               ##
#                       By default this file is called 'ParamFiles.txt' and   ##
#                       -if the full path it is not specified- it is searched ##
#                       for within the 'PSO.in' subdirectory of 'drty.model'  ## 
# 'param.ranges'      : character with the name of the hydroPSO input file    ##
#                       defining the minimum and maximum boundary values for  ##
#                       each one of the parameters to be calibrated.          ##
#                       By default this file is called 'ParamRanges.txt' and  ##
#                       -if the full path it is not specified- it is searched ##
#                       for within the 'PSO.in' subdirectory of 'drty.model'  ## 
# 'observations.fname': character with the name of the hydroPSO output file   ##
#                       storing the observed values used during the optimisa- ##
#                       tion.                                                 ##
#                       By default this file is called 'Observations.txt' and ##
#                       -if the full path it is not specified- it is searched ##
#                       for within the  'PSO.out'                             ##
#                       subdirectory of 'drty.model'.                         ##
# 'exe.fname'         : character, model command line arguments to be entered ##
#                       through a prompted string to execute the user-defined ##
#                       model                                                 ## 
# 'decimals'          :                                                       ##
################################################################################
pest2hydroPSO <- function(pst.fname, 
                          drty.pest=NULL, 
                          drty.model=NULL, 
                          drty.out="PSO.in",
                          param.files="ParamFiles.txt",
                          param.ranges="ParamRanges.txt",
                          decimals=5,
                          verbose=TRUE) {
   
  # Checkings
  if (missing(pst.fname)) stop("PEST control file is missing ('pst.fname')")                      
  if (is.null(drty.pest)) drty.pest <- dirname(pst.fname)
  if (is.null(drty.model)) drty.model <- dirname(pst.fname)
  
  if (drty.out=="PSO.in") drty.out <- paste(dirname(pst.fname), "/PSO.in", sep="")
  if (!file.exists(drty.out)) dir.create(drty.out, recursive=TRUE)
  
  if (basename(param.files)==param.files) 
    param.files <- paste(drty.out, "/", param.files, sep="")
    
  if (basename(param.ranges)==param.ranges) 
    param.ranges <- paste(drty.out, "/", param.ranges, sep="")

  ##############################################################################
  # 1) Reading .pst file
  if (verbose) message("                                                  ")
  if (verbose) message("[ 1) Reading the .pst file '", pst.fname, "' ... ]")
  x <- readLines(pst.fname)
  
  ##############################################################################
  # Getting the number of parameters and observations
  values    <- as.numeric(strsplit(x[4], " ")[[1]])
  nna.index <- which(!is.na(values))
  
  nparam <- values[nna.index][1]
  nobs   <- values[nna.index][2]

  if (verbose) message("                                                 ")
  if (verbose) message("   [ Number of parameters found  : ", nparam, " ]")
  if (verbose) message("   [ Number of observations found: ", nobs, " ]")  
  
  # Getting the number of input files (.tpl) and output files (.ins)
  files     <- strsplit(x[5], " ")[[1]]
  spc.index <- which(files=="")  
  if (length(spc.index) > 0) files <- files[-spc.index]
  ntpl      <- as.numeric(files[1])
  nins      <- as.numeric(files[2])

  if (verbose) message("   [ Number of .tpl found        : ", ntpl, " ]") 
  if (verbose) message("   [ Number of .ins found        : ", nins, " ]")    
  
  # Getting Param names and Ranges
  L <- 0
  params.stg <- "* parameter data"
  ini.index <- which(x==params.stg) 
  
  L <- length(ini.index)
  if (L > 0) {
    suppressWarnings(tmp <- read.table(file=pst.fname,skip=ini.index,header=FALSE,nrows=nparam,colClasses
        =c("character","NULL","NULL","numeric","numeric","numeric","NULL","NULL","NULL","NULL"),fill=TRUE,na.strings=""))
        
    param.names <- tmp[,1]
    param.ini   <- as.numeric(tmp[,2])
    param.min   <- as.numeric(tmp[,3])
    param.max   <- as.numeric(tmp[,4])
  } else stop("Invalid pst file: ", params.stg, " does not exist !")
  
  ##############################################################################
  # 2) Writing 'ParamRanges.txt'
  if (verbose) message("                                                          ")
  if (verbose) message("[ 2) Writing 'ParamRanges.txt' into '", drty.out, "' ... ]")
  .pst2paramranges(drty.model=drty.model, names=param.names, ini=param.ini,
                  min=param.min, max=param.max, 
                  fname.out=param.ranges)
  
  
  ##############################################################################
  # 3) Writing file with observations
  if (verbose) message("                                                                ")
  if (verbose) message("[ 3) Writing 'PEST2hydroPSO_OBS.txt' into '", drty.out, "' ... ]")
  L <- 0
  obs.stg <- "* observation data"
  ini.index <- which(x==obs.stg) 
  
  L <- length(ini.index)
  if (L > 0) {
    suppressWarnings(tmp <- read.table(file=pst.fname,skip=ini.index,header=FALSE,nrows=nobs,colClasses
        =c("NULL","numeric","NULL","NULL"),fill=TRUE,na.strings=""))
        
    obs <- tmp[,1]
  } else stop("Invalid pst file: ", obs.stg, " does not exist !")
  
  # Writing the output file
  obs.fname <- paste(drty.out, "/PEST2hydroPSO_OBS.txt", sep="")
  write.table(obs, file=obs.fname, col.names=FALSE, row.names=FALSE, quote=FALSE)
  

  ##############################################################################  
  # 4) Model command line
  if (verbose) message("                                       ")
  if (verbose) message("[ 4) Searching model command line ... ]")
  L <- 0
  model.stg <- "* model command line"
  ini.index <- which(x==model.stg) 
  
  L <- length(ini.index)
  if (L > 0) {        
    model.exe <- x[ini.index+1]
  } else stop("Invalid pst file: ", model.stg, " does not exist !")

  if (verbose) message("   [ Model command line          : '", model.exe, "' ... ]")
  

  ##############################################################################  
  # 5) Getting Parameter filenames and locations (.tpl's and .ins's)
  if (verbose) message("                                     ")
  if (verbose) message("[ 5) Writing model input/output ... ]")
  if (verbose) message("                                     ")
  
  L <- 0
  io.stg <- "* model input/output"
  ini.index <- which(x==io.stg) 
  
  L <- length(ini.index)
  if (L > 0) {
    suppressWarnings(tmp <- read.table(file=pst.fname,skip=ini.index,header=FALSE,nrows=ntpl+nins,colClasses
        =c("character","character"),fill=TRUE,na.strings=""))
        
    tpls   <- tmp[1:ntpl,1]
    inputs <- tmp[1:ntpl,2]
    ins    <- tmp[(ntpl+1):(ntpl+nins),]
  } else stop("Invalid pst file: ", io.stg, " does not exist !")
  
  # Writing ParamFiles.txt
  .pst2paramfiles(drty.model=drty.model, tpls=tpls, inputs=inputs, 
                 param.names=param.names, fname.out=param.files, 
                 decimals=decimals)
                 
                 
  ##############################################################################  
  # 6) Creating the Rscript used for running hydroPSO
  if (verbose) message("                                                ")
  if (verbose) message("[ 6) Copying R script for running hydroPSO ... ]")
  rscript.fname <- system.file("Rscripts/hydroPSO-Rscript.R", package="hydroPSO")
  dst.fname     <- paste(drty.model, "/Rscripts/hydroPSO-Rscript.R", sep="")
  if (!file.exists(dirname(dst.fname))) dir.create(dirname(dst.fname), recursive=TRUE)
  file.copy(rscript.fname, dst.fname, overwrite=TRUE)
  
  ##############################################################################  
  # 7) Modifying the Rscript used for running hydroPSO
  if (verbose) message("                                                  ")
  if (verbose) message("[ 7) Modifying R script for running hydroPSO ... ]")
  
  # rading the Rscript
  x <- readLines(dst.fname)
  
  # Changing model directory
  x[25] <- paste("model.drty <- \"", drty.model, "\"", sep="")
  
  # Changing model exe name
  x[64] <- paste("exe.fname= \"", model.exe, "\"", sep="")  
  
  # Writing the modified file
  if (file.exists(dst.fname)) file.remove(dst.fname)
  writeLines(x, con=dst.fname)
  
  ##############################################################################  
  # 8) Output
  if (verbose) message("                                                           ")
  if (verbose) message("[          PEST2hydroPSO finished !!                      ]")
  if (verbose) message("[ R script to run hydroPSO available in: '", dst.fname, "']")
  if (verbose) message("[ Before running hydroPSO, you MUST modify the section:    ")
  if (verbose) message("  'User-defined variables' in '", basename(dst.fname), "'")
  
} # 'pest2hydroPSO' END
