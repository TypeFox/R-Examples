# Part of the hydroPSO package, http://www.rforge.net/hydroPSO/
# Copyright 2010-2014 Mauricio Zambrano-Bigiarini
# Distributed under GPL 2 or later

################################################################################
#                             ModifyInputfile                                  #
################################################################################
# Purpose    : To write a numeric value into a specified position of a plain   #
#              text file                                                       #
################################################################################
# Output     : A mofified text file ('filename')                               #
################################################################################
# Author : Mauricio Zambrano-Bigiarini                                         #
# Started: 17-Dec-2010 at JRC Ispra                                            #
# Updates: 20-Jan-2011                                                         #
#          06-Sep-2013                                                         #
#          09-Abr-2014                                                         #
################################################################################

ModifyInputFile <- function(
                            ParamID,    # character, with the ID of the parameter to be modified (only used for visualization purposes)
                            newvalue,   # numeric value to be written into the text file
                            filename,   # character, with the name of the text file that will be modified
                            row,        # numeric, with the row number in \code{filename} where \code{newvalue} will be written
                            col.ini,    # numeric, with the starting column number in \code{filename} where \code{newvalue} is going to be written.
                            col.fin,    # numeric, with the ending column number in \code{filename} where \code{newvalue} is going to be written.
                            decimals,   # numeric, with the number of decimal places used to write \code{newvalue} into \code{filename}
                            verbose=TRUE) {

  if (!file.exists(filename))
    stop( paste("Invalid argument: the file '", filename, "' doesn't exist!", sep="") )

  lines  <- readLines(filename)

  myline <- lines[row]

  L.trg <- col.fin - col.ini + 1

  newvalue.stg <- as.character(round(newvalue, decimals))

  L <- nchar(newvalue.stg)
  
#  message("ParamID  : ", ParamID)
#  message("filename : ", basename(filename))
#  message("row      : ", row)
#  message("new value: ", newvalue.stg)
#  message("L.trg    : ", L.trg)

  if (L < L.trg) newvalue.stg <- format(newvalue, justify="right", width=L.trg, nsmall=decimals)  
  
  if (L > L.trg) {
     #newvalue.stg <- format(newvalue, justify="right", width=L.trg, scientific=TRUE)
     #e.pos <- which(strsplit(newvalue.stg, split=character(0))[[1]] == "e")
     #newvalue.stg <- format(newvalue, justify="right", width=L.trg, scientific=TRUE, digits=e.pos-1)
     nexp <- 2
     if (abs(newvalue) >= 1E100) nexp <- 3
     dig          <- max(decimals-(L - L.trg)-3-nexp, 0) 
     suppressWarnings(
       newvalue.stg <- formatC(newvalue, width=L.trg, format="E", digits=dig)
     )
  } # IF end 
   
  substr(myline, col.ini, col.fin) <- newvalue.stg

  lines[row] <- myline

  writeLines(lines, filename)

  if (verbose)
   message( paste("[", ParamID, ": '", round(newvalue,5), "' was successfully put into '", basename(filename), "']", sep="") )

} # 'ModifyInputFile' END
