################################################################################
#                           read.paramfile                                     #
################################################################################
# Purpose    : Reading 'ParamFiles.txt' and 'ParamRanges.txt' files            #
################################################################################
# Output     : A data.frame with the following columns                         #
#              For 'ParamRanges.txt' file:                                     #
#              1) ParameterNmbr: numeric, from 1 to the number of parameters   #
#              2) ParameterName: character, with a meaningful name of each parameter
#              3) MinValue     : numeirc, with the minum value that each parameter can take 
#              3) MaxValue     : numeric, with the maximum value that each parameter can take 
#
#              For 'ParamFiles.txt' file:                                      #
#              1) ParameterNmbr: numeric, from 1 to the number of parameters   #
#              2) ParameterName: character, with a meaningful name of each parameter
#              3) Filename     : character, with the name of each file that has to be modified
#              4) Row.Number   : numeric, with the number of the row in 'ParamFiles.txt' that has to be modified
#              5) Col.Start    : numeric, with the number of the colum in 'ParamFiles.txt' where the modification starts
#              6) Col.End      : numeric, with the number of the colum in 'ParamFiles.txt' where the modification ends
#              7) Decimals     : numeric, with the number of decimal places to be used for writing the new values in 'ParamFiles.txt'
################################################################################
# Author     : Mauricio Zambrano-Bigiarini                                     #
# Started    : 13-Dec-2010 at JRC Ispra                                        #
# Last Update: 13-Dec-2010                                                     #
################################################################################
# Function for 
read.paramfile <- function(file="ParamFiles.txt") {

  x <- read.table(file=file, header=TRUE, skip=0, stringsAsFactors=FALSE)
 
  return(x)
 
} #'read.paramfile' END
