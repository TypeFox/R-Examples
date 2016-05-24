trimToTaxaLevel <-
function(data, level="genus", eliminateParentNodes=FALSE, trimBelow=NULL, split="."){
if(missing(data))
stop("A valid data set is required.")

maxLevel <- getTaxaDepth(level)

nameSplit <- strsplit(as.character(rownames(data)), split, fixed=TRUE)
lowerLevels <- NULL
for(l in 1:nrow(data)){ 
if(length(nameSplit[[l]]) == maxLevel){
lowerLevels <- c(lowerLevels, l)
if(eliminateParentNodes && is.null(trimBelow))
rownames(data)[l] <- nameSplit[[l]][maxLevel]
}else if(!is.null(trimBelow) && is.logical(trimBelow)){
if(length(nameSplit[[l]]) < maxLevel && trimBelow){
lowerLevels <- c(lowerLevels, l)
}else if(length(nameSplit[[l]]) > maxLevel && !trimBelow){
lowerLevels <- c(lowerLevels, l)
}
}
}
data <- data[lowerLevels, , drop=FALSE]

return(data)
}
