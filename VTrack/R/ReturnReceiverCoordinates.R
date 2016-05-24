ReturnReceiverCoordinates <-
function(sReceiver, iArraySize, CoordArray)
{
  sResult <- ""
  for (i in 1:iArraySize)
    if (sReceiver == CoordArray[i,1])
      # flip the lat/long to long/lat for KML
      sResult <- paste(CoordArray[i,3],CoordArray[i,2],sep=",")
  
  sResult
}
