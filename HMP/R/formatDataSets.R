formatDataSets <-
function(data){
if(missing(data))
stop("data missing.")

numData <- length(data)
if(numData < 2)
stop("At least 2 data sets are required.")

newData <- NULL
for(i in 1:length(data)){
newData <- merge(newData, t(data[[i]]), by=0, all=TRUE)
rownames(newData) <- newData[,1]
newData <- newData[,-1]
}
newData[is.na(newData)] <- 0
newData <- t(newData)
newData <- newData[,colSums(newData) != 0]
newData <- newData[rowSums(newData) != 0,]
newData <- newData[,order(colSums(newData), decreasing=TRUE)]

retData <- vector("list", numData)
base <- 0
for(i in 1:numData){
retData[[i]] <- newData[(base+1):(nrow(data[[i]])+ base),]
base <- base + nrow(data[[i]])
}

names(retData) <- names(data)
return(retData)
}
