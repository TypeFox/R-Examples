# Part of the hydroPSO package, http://www.rforge.net/hydroPSO/
# Copyright 2010-2012 Mauricio Zambrano-Bigiarini
# Distributed under GPL 2 or later

################################################################################
#                           read.ParameterRanges                               #
################################################################################
# Purpose: To read a file containing the minimum and maximum values within     #
#          one or more parameter can be varied within an optimisation          #
#          algorithm                                                           #
################################################################################
# Output : A data.frame with the following columns                             #
#          1) ParameterNmbr: numeric, from 1 to the number of prameters        #
#          2) ParameterName: character, with a meaningful name of each parameter
#          3) MinValue     : numeric, with the minum value that each parameter can take 
#          4) MaxValue     : numeric, with the maximum value that each parameter can take 
################################################################################
# Author     : Mauricio Zambrano-Bigiarini                                     #
# Started    : 13-Dec-2010 at JRC Ispra                                        #
# Last Update: 13-Dec-2010 ; 11-Jan-2012                                       #
################################################################################
read.ParameterRanges <- function(ParamRanges.fname="ParamRanges.txt" # It has to have a row with the header
                                 ) {

  # Reading the file with the Range of the paramters
  ParamRanges <-  read.paramfile(file=ParamRanges.fname)

  n1 <- length(levels(as.factor(ParamRanges[,1]))) # Number of Param IDs
  n2 <- length(levels(as.factor(ParamRanges[,2]))) # Number of Param Names
  if ( n1 != n2)
    stop( paste("In '", basename(ParamRanges.fname), "' : Number of 'ParameterNmbr' != Number of 'ParameterName' (", n1, " != ", n2, ")", sep="" ) ) 
  
  # Number of parameters that have to be optimised
  nparam <- nrow(ParamRanges)

  message( paste("[ Number of parameters:", nparam, " ]", sep=" ") )
  message( paste("[ Parameters' names   : ", paste(ParamRanges[,2], collapse = ", "), " ]", sep="" ) )

  # Giving a meaningful name to each parameter in 'X.MinMax'
  rownames(ParamRanges) <- ParamRanges[,2]

  # Returning only the 3rd and 4th column, with the minim and maximum vpossible values, respectively
  return(ParamRanges[,3:4])

} # 'read.ParameterRanges' END

#read.ParameterRanges("./PSO/ParamRanges.txt")
