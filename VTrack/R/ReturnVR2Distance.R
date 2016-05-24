ReturnVR2Distance <-
function(NonResidenceFile,sDistanceMatrix)
{
  iReceiverIndex1 <- 0
  iReceiverIndex2 <- 0
  VR2Distance <- 0
  
  for(j in 1:dim(NonResidenceFile)[1])
  {
    for (i in 1:dim(sDistanceMatrix)[1])
    {    
      if (as.character(sDistanceMatrix$DM[i]) == as.character(NonResidenceFile$RECEIVERID1[j]))
        iReceiverIndex1 <- i
      if (as.character(sDistanceMatrix$DM[i]) == as.character(NonResidenceFile$RECEIVERID2[j]))
        iReceiverIndex2 <- i
    }
    VR2Distance[j] <- sDistanceMatrix[iReceiverIndex1,(iReceiverIndex2)+1]
  }
  
  return(VR2Distance)
}
