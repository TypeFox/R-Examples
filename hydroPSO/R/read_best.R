# File read_best.R
# Part of the hydroPSO R package, http://www.rforge.net/hydroPSO/ ; 
#                                 http://cran.r-project.org/web/packages/hydroPSO
# Copyright 2011-2012 Mauricio Zambrano-Bigiarini & Rodrigo Rojas
# Distributed under GPL 2 or later

################################################################################
#                            'read_best'                                       # 
################################################################################
# Purpose:                                                                     #
# To read the 'BestParameterSet.txt' ouput file, wich stores the best parameter#
# set and its corresponding best fitness value                                 #
################################################################################
# Author : Mauricio Zambrano-Bigiarini & Rodrigo Rojas                         #  
# Started: 08-Nov-2011,                                                        #
# Updates:                                                                     #        
################################################################################

read_best <- function(file="BestParameterSet.txt",
                      verbose=TRUE
                      ) {
                         
                         
  # Checking that 'file' exists
  if ( !file.exists(file) )
     stop( paste("Invalid argument value: The file '", basename(file), "' doesn't exist", sep="" ) )

  # Reading the file
  if (verbose) message( "                                                     ")  
  if (verbose) message( "[ Reading the file '", basename(file), "' ... ]" )  
  best  <- read.table(file=file, header=TRUE, skip=0)

  # computing the number of columns
  ncols <- ncol(best)  
  
  # Getting the number of the 'best particle'
  best.part.number <- best[, 1]

  # Getting the GoF of the 'best particle'
  best.param.gof <- best[, 2]
  
  # Getting the Parameter Values of the 'best particle'
  best.param.values <- best[, 3:ncols]
  
  out <- list(best.part.number=best.part.number, 
              best.param.values=best.param.values, 
              best.param.gof=best.param.gof)
  
  return(out)
  
}  # 'read_best' END
