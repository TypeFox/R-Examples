################################################################
#' Open file to write result
#' @param filename file name
#' Output: Results in text file 
#' @export
################################################################
# DECODE (c) Copyright 2014 by The Hong Kong Polytechnic University, Department of Health Technology and Informatics 
# Written by Thomas Lui
# Permission is granted to copy and use this program provided no fee is
# charged for it and provided that this copyright notice is not removed.
##########################################################################
openFileToWrite = function(filename) {

  outFile = file(filename,"w")
  
  #########################################
  # printing result to the database
  # To simplify the results for basic user, some detail results are commented
  #########################################

  writeLines("id", outFile,sep="\t")
  writeLines("Gene", outFile,sep="\t")
  writeLines("t score", outFile,sep="\t")
  writeLines("abs t score", outFile,sep="\t")
  writeLines("p-value", outFile,sep="\t")
  writeLines("Fold change", outFile,sep="\t")

  writeLines("optimal_Tcutoff", outFile,sep="\t")
  writeLines("optimal_DRcutoff", outFile,sep="\t")
  writeLines("r+,t+", outFile,sep="\t")
  writeLines("r+,t-", outFile,sep="\t")
  writeLines("r-,t+", outFile,sep="\t")
  writeLines("r-,t-", outFile,sep="\t")
  writeLines("expect(r+,t+)", outFile,sep="\t")
  writeLines("chi square value", outFile,sep="\t")
  writeLines("p-value", outFile,sep="\t")

  writeLines("Best associated gene set", outFile,sep="\t")
  writeLines("category", outFile,sep="\t")
  writeLines("tempA", outFile,sep="\t")
  writeLines("expected A", outFile,sep="\t")
  writeLines("p-value", outFile,sep="\t")
  writeLines("matched genes", outFile,sep="\t")
  writeLines("FR in low DE", outFile,sep="\t")
  writeLines("FR in high DE", outFile,sep="\t")
  writeLines("FR in low r", outFile,sep="\t")
  writeLines("FR in high r", outFile,sep="\t")
  writeLines("FR in low r, low DE", outFile,sep="\t")
  writeLines("FR in low r, high DE", outFile,sep="\t")
  writeLines("FR in high r, low DE", outFile,sep="\t")
  writeLines("FR in high r, high DE", outFile,sep="\t")
 

  writeLines("",outFile,sep="\n")
  close(outFile)
}