formatData <-
function(data, countThreshold=1000, normalizeThreshold=10000){
if(missing(data))
stop("A valid data set is required.")

data <- data[order(rownames(data)),]
data[is.na(data)] <- 0
data <- data[, data[1,] >= countThreshold]

if(ncol(data) == 0)
stop("'countThreshold' is too high.")

if(normalizeThreshold > 0){
for(i in ncol(data):1)
data[,i] <- data[,i] * (normalizeThreshold/data[1,i])
}

return(data)
}
