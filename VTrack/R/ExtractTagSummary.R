ExtractTagSummary <-
function(sInputFile)
{
  TransmitterList <- ExtractUniqueValues(sInputFile,2)

  fExtractsummaryid <- function(i){
    T1 <- ExtractData(sInputFile,sQueryTransmitterList = TransmitterList[i])
    row.names(T1)<- NULL
    TRANSMITTERID <- T1[1,2]
    STARTDATE <- T1[1,1]
    STARTREC <- T1[1,5]
    NODETECTS <- nrow(T1)
    ENDDATE <- T1[NODETECTS,1]
    ENDREC <- T1[NODETECTS,5]
    NORECS <- length(unique(T1$RECEIVERID))
    return(data.frame(TRANSMITTERID,STARTDATE,ENDDATE,NODETECTS,STARTREC,ENDREC,NORECS))
  }
  
  lnewdf <- lapply(1:length(TransmitterList),fExtractsummaryid)
  newdf <- do.call(rbind,lnewdf)
  return(newdf)
}
