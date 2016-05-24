Data.filter <-
function(data, order.type, reads.crit, K){
if(missing(data) || missing(order.type) || missing(reads.crit) || missing(K))
stop("data, order.type, reads.crit and/or K missing.")

if(K >= ncol(data))
stop(sprintf("K is too large.  It must be smaller than %i.", ncol(data)-1))
if(K <= 0)
stop("K is too small.  It must be larger than 0.")
if(reads.crit >= max(rowSums(data)))
stop(sprintf("Reads.crit is too large.  It must be smaller than %i.", max(rowSums(data))))

data <- data[apply(data, 1, sum)>reads.crit, ,drop=FALSE]

if(tolower(order.type) == "sample"){
data <- t(apply(data, 1, function(x){x[order(x, decreasing=TRUE)]}))
}else if(tolower(order.type) == "data"){
data <- data[,order(colSums(data), decreasing=TRUE)]
}else{
data <- data[,order(colSums(data), decreasing=TRUE)]
warning(sprintf("order.type defaulting to 'data'. '%s' not recognized.", as.character(order.type)))
}

if(nrow(data) < 2)
stop("Reads.crit is so large that it excludes all but one sample.  Please try lowering its value.")

data <- data[,apply(data, 2, sum)>0]
data <- cbind(data[,1:K], as.matrix(apply(as.matrix(data[,(K+1):ncol(data)]), 1, sum)))
colnames(data)[K+1] <- "Others"

return(data)
}
