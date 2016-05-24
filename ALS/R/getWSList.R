`getWSList` <-
function(S, WList, tt) {
  SList <- vector("list", length=length(WList))
  for(j in 1:length(WList)) 
    SList[[j]] <-  t(S * WList[[j]][tt,])
  SList 
}

