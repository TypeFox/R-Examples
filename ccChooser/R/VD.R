VD <-
function(CC, EC){
var_CC <- apply(CC, 2, var)
var_EC <- apply(EC, 2, var)
tmp <- unlist(lapply(1:ncol(CC), function(i, var_CC, var_EC){
return( ( var_CC[i] - var_EC[i] ) / var_EC[i] )
}, var_CC, var_EC))
VD_index <- ( sum(tmp) / ncol(CC) ) * 100

return(VD_index)
}
