LOG <-
function(groups, size, newSize, grpIDs, grpSize){
sumGrpLog<-unlist( lapply( 1:length(grpIDs), function(i, grpSize){
return(log(grpSize[i]) * grpSize[i])
}, grpSize) )
sumLog <- sum(sumGrpLog)

newGrpSizes<-unlist( lapply( 1: length(grpIDs), function(i, newSize, sumGrpLog, sumLog){
newGrpSize <- round(newSize * ( sumGrpLog[i] / sumLog ))
if(newGrpSize == 0)
newGrpSize <- 1
return(newGrpSize)
}, newSize, sumGrpLog, sumLog) )

return(newGrpSizes)
}
