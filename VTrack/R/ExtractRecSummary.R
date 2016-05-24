ExtractRecSummary <-
function(sInputFile)
{
  ReceiverList <- ExtractUniqueValues(sInputFile,5)
  
  fExtractsummaryRecid <- function(i){
    T1 <- ExtractData(sInputFile,sQueryReceiverList = ReceiverList[i])
    T1 <- T1[order(T1$DATETIME),]
    row.names(T1)<- NULL
    RECEIVERID <- T1[1,5]
    STATIONNAME <- T1[1,6]
    FIRSTDETECT <- min(T1[1,1])
    NODETECTS <- nrow(T1)
    LASTDETECT <- T1[NODETECTS,1]
    NOTRANSMITTER <- length(unique(T1$TRANSMITTERID))
    return(data.frame(RECEIVERID,STATIONNAME,FIRSTDETECT,NODETECTS,LASTDETECT,NOTRANSMITTER))
  }

  lnewdf <- lapply(1:length(ReceiverList),fExtractsummaryRecid)
  newdf <- do.call(rbind,lnewdf)
  return(newdf)
}
