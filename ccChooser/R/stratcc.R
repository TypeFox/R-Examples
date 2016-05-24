stratcc <-
function(x, groups, alloc = "Pro", fraction = 0.1, clustering = FALSE, cluster_method = "ward"){
allocated <- allocc(x, groups, fraction, method = alloc)
grpIDs <- unique(groups)
if(clustering){
CC <- lapply(grpIDs, function(grpID){
grpData <- subset(x, groups == grpID)
if(nrow(grpData) == 1){
return(grpData)
} else {
tree <- hclust( daisy(grpData, stand = TRUE), method = cluster_method)
clGrpID <- cutree(tree, k = allocated[ grpID, 2 ])
tmp <- lapply(1:allocated[allocated[,1] == grpID, 2], function(alloID){
alloData <- grpData[clGrpID == alloID , ]
resID <- sample(1:nrow(alloData), 1)
return(alloData[resID , ])
})
}
tmp <- as.data.frame(do.call(rbind, tmp), stringsAsFactors = FALSE)
return(tmp)
})
}else{
CC <- lapply(grpIDs, function(grpID){
grpData <- subset(x, groups == grpID)
allocID <- sample( 1:nrow(grpData), allocated[ allocated[,1]==grpID, 2 ] )
newData <- grpData[ allocID, ]
return(newData)
})
}
CC<-as.data.frame(do.call(rbind, CC), stringsAsFactors = FALSE)
return(CC)
}
