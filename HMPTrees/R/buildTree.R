buildTree <-
function(data, level="genus", split="."){
if(missing(data))
stop("A valid data set is required.")

names <- rownames(data)
nameSplit <- strsplit(as.character(names), split,  fixed=TRUE)

retData <- data

for(lvl in 1:(getTaxaDepth(level)-1)){
newData <- data.frame()
for(i in 1:length(nameSplit)){ 
name <- paste(nameSplit[[i]][1:lvl], collapse=split)
if(is.element(name, rownames(newData))){
loc <- grep(name, rownames(newData), fixed=TRUE)
newData[loc,] <- newData[loc,] + data[i,]
}else{
newData <- rbind(data[i,], newData)
rownames(newData)[1] <- name
}
}
retData <- rbind(retData, newData)
}

retData <- ifelse(retData>=1, 1, 0)
retData <- retData[order(rownames(retData)),]

return(retData)
}
