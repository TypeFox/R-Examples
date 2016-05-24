removeUnclass <-
function(data, remove=TRUE){
falsePos <- grep(".U", rownames(data), fixed=TRUE)
unclass <- grep("U", rownames(data), fixed=TRUE)

if(length(unclass) == 0) #nothing to remove
return(data)

if(length(falsePos) != 0) #we have false positives
unclass <- setdiff(unclass, falsePos)

if(remove){
data <- data[-unclass,]
}else{
data[unclass,] <- 0
}

return(data)
}
