# This script file takes an annotation data frame (Format: Name, Chr, Start, Stop) and then
# a dataframe that specifies regions of interest. THe output is then just the subset of
# this

getRegionsOI <- function(annot,regOI){
  tempResult <- data.frame(Name="-1",Chr="-1",Start=-1,End=-1)
  for(oiRun in 1:nrow(regOI)){
    temp <- annot[annot[,2]==regOI[oiRun,1],]
    colnames(temp) <- c("Name","Chr","Start","End")
    temp <- temp[temp[,3]>=regOI[oiRun,2],]
    temp <- temp[temp[,3]<=regOI[oiRun,3],]
    tempResult <- rbind(tempResult,temp)
  }
  tempResult <- tempResult[-1,]
  tempResult[,1] <- as.character(tempResult[,1])
  tempResult
}