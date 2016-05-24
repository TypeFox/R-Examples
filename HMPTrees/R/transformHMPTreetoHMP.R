transformHMPTreetoHMP <-
function(data, elimZero=FALSE, zeroValue=.00001){
if(missing(data))
stop("A valid data set is required.")

dataSum <- rowSums(data)
loc <- NULL
for(r in 1:nrow(data)){
if(dataSum[r] == 0){
loc <- c(loc, r)
data[r, 2] <- zeroValue
}
}
if(elimZero)
data <- data[-loc,]

return(t(data))
}
