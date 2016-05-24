D3 <-
function(size, newSize, grpIDs, grpSize, groups, data){
sumGrpDist <- unlist( lapply( 1: length(grpIDs), function(i){
tmp <- mean( daisy( data[groups==grpIDs[i] , ] ) )* log(grpSize[i]) * grpSize[i]
if(is.nan(tmp))
tmp <- 0
return(tmp)
}))

sumDist <- sum(sumGrpDist)

newGrpSizes <- unlist( lapply( 1: length(grpIDs), function(i){
newGrpSize <- round(newSize * (sumGrpDist[i] / sumDist))
if(newGrpSize == 0)
newGrpSize <- 1
return(newGrpSize)
}) )

return(newGrpSizes)
}
