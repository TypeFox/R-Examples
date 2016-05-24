######################################################################
#
# toSPIAData: convert SNPs data into SPIA format
#
# Input: 1. SNPMatrix: a matrix with a column for each cell line and a 
#          row for each SNP
#        2. encoding: a four elements ecnoding vector describing the 
#          encoding used in SNPMatrix. 
#          For instance, (0,2,1,-1) says that SNPMatrix uses 0 for AA,
#          2 for BB, 1 for AB, and -1 for NoCall
# Output: 1. a matrix with a column for each cell line and a 
#           row for each SNP encoded with SPIA encoding
#
######################################################################

toSPIAData <- function(SNPMatrix, encoding){
  
  ##SPIA encoding vector
  SPIA_values <- as.integer(c(0,1,2,NA))
  
  ##Four constants vector for temporary substitutions
  SPIA_tmp <- c("SPIA_tmp_AA", "SPIA_tmp_BB", "SPIA_tmp_AB","SPIA_tmp_NoCall")
  
  
  ##Control the lenght of the encoding vector
  if (!length(encoding) == 4){
    message("SPIA: error in function toSPIAData encoding vector must contain four elements.");
    return(-1)
  }  
  
  ##A temporary substitution is needed to avoid clash of values
  
  ##Replace all the value of SNPMatrix with temporary constants
  for( i in c(1:4) ) SNPMatrix <- replace(SNPMatrix, SNPMatrix == encoding[i], SPIA_tmp[i])
  
  ##Replace temporary constansts with SPIA values
  for( i in c(1:4) ) SNPMatrix <- replace(SNPMatrix, SNPMatrix == SPIA_tmp[i], SPIA_values[i])
  
  message("SPIA: encoding completed.");
  return(SNPMatrix)
  
}