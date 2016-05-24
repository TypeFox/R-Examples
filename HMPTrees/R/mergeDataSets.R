mergeDataSets <-
function(data, calcMLE=TRUE, uniqueNames=FALSE){
newData <- NULL
mles <- NULL

if(missing(data))
stop("At least 1 valid data set is required.")

if(class(data) == "data.frame" || class(data) == "matrix"){ #turn a single dataset into a list
temp <- data
data <- NULL
data[[1]] <- temp
}

for(i in 1:length(data)){
if(uniqueNames)
colnames(data[[i]]) <- paste("Data", i, "-", colnames(data[[i]]))
newData <- merge(newData, data[[i]], by=0, all=TRUE)
rownames(newData) <- newData[,1]
newData <- newData[,-1]

if(calcMLE)
mles[[i]] <- getMLEandLoglike(data[[i]])$mleTree
}
newData[is.na(newData)] <- 0 #set any NA's to 0

if(calcMLE){
mleAll <- getMLEandLoglike(newData)$mleTree
colnames(mleAll) <- "combined mle"

if(length(data) > 1){
for(i in 1:length(mles)){
newData <- merge(newData, mles[[i]], by=0, all=TRUE)
rownames(newData) <- newData[,1]
newData <- newData[,-1]

colnames(newData)[ncol(newData)] <- paste("data", i, "mle", sep="")
}
newData[is.na(newData)] <- 0 #set any NA's to 0
}

newData <- cbind(newData, mleAll)
}

return(newData)
}
