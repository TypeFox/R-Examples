ConvertPointsToCircuitousDistance <-
function(sInputFile)
{
  iCount <- 0
  iLog <- 0
  rDistance <- 0
  rPreviousLat <- 0
  rPreviousLon <- 0
  sPreviousRECEIVERID <- ""
  
  # create output files 
  outfile <- data.frame(RECEIVERID1=sInputFile[1,1],RECEIVERID2=sInputFile[1,1],DISTANCE=0)[NULL,]
  
  while (iCount < nrow(sInputFile))
  {
    # increment line counter
    iCount <- iCount + 1
    sRECEIVERID <- sInputFile[iCount,1]
    rLat <- sInputFile[iCount,2]
    rLon <- sInputFile[iCount,3]
    
    if (iCount > 1)
    {
      # compute distance from previous point
      rSegmentDistance <- ComputeDistance(rLat,rPreviousLat,rLon,rPreviousLon)
      rDistance <- rDistance + rSegmentDistance
      if (sRECEIVERID != 0)
      {
        iLog <- iLog + 1 
        # write record to outfile
        outfile[iLog,] <- c(sPreviousRECEIVERID,sRECEIVERID,rDistance)
        rDistance <- 0
      }
    }
    rPreviousLat <- rLat
    rPreviousLon <- rLon
    if (sRECEIVERID != 0)
    {
      sPreviousRECEIVERID <- sRECEIVERID
    }
  }
  return(outfile)
}
