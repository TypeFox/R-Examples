RR <-
function(CC, EC){
colDiff_CC <- apply(CC, 2, max) - apply(CC, 2, min)
colDiff_EC <- apply(EC, 2, max) - apply(EC, 2, min)
tmp <- unlist(lapply(1:ncol(CC), function(i, colDiff_CC, colDiff_EC){
return( (colDiff_CC[i] - colDiff_EC[i]) / colDiff_EC[i] )
}, colDiff_CC, colDiff_EC))
RR_index <- (sum(tmp) / ncol(CC)) * 100

return(RR_index)
}
