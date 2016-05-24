`getWSListAllTimes` <-
function(S,WList) {
  allWSList <- vector("list",nrow(WList[[1]]))
  for(i in 1:nrow(WList[[1]])) 
    allWSList[[i]] <- getWSList(S, WList, i)
  allWSList
}

