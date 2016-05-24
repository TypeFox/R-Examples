allocc <-
function(x, groups, fraction=0.1, method="Pro"){
if(!any(method %in% c("Pro", "Log", "D2", "D3")))
stop("Bad name for allocation method!")
size <- nrow(x)
if(size!=length(groups))
stop("Number of groups must be equal number of rows in x")
newSize <- size*fraction

grpIDs <- sort(unique(groups))
grpSize<-unlist(lapply(grpIDs, function(i, groups){sum(groups==i)}, groups))


if (method=="Pro"){

cat("Proportional allocation method\n")
newGrpSizes<-PRO(groups, size, newSize, grpIDs)

} else if (method=="Log"){

cat("Logarytmic allocation method\n")
newGrpSizes<-LOG(groups, size, newSize, grpIDs, grpSize)

} else if (method=="D2"){

cat("D2 allocation method\n")
newGrpSizes <- D2(size, newSize, grpIDs, grpSize, groups, x)

} else if (method=="D3"){

cat("D3 allocation method\n")
newGrpSizes <- D3(size, newSize, grpIDs, grpSize, groups, x)

}

results<-cbind(grpIDs,newGrpSizes, deparse.level=0)
cat(paste("Number of accessions in core colection:", sum(results[,2]),"\n" ,sep=" "))

results <- results[ order(results[ , 1]), ]
colnames(results) <- c("GroupID", "NewSize")
invisible(results)
}
