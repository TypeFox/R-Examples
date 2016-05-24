ReturnReceiverIndex <-
function(sReceiver, iArraySize, CoordArray)
{
  iResult <- -1
  for (i in 1:iArraySize)
    if (sReceiver == CoordArray[i,1])
      iResult <- i
  
  iResult
}
