PRO <-
function(groups, size, newSize, grpIDs){
newGrpSizes<-unlist(lapply(1:length(grpIDs), function(i, groups){
newGrpSize <- round(newSize*(sum(groups==i)/size))
if(newGrpSize == 0)
newGrpSize <- 1
return(newGrpSize)
}, groups))

return(newGrpSizes)
}
