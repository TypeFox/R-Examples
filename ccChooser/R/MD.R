MD <-
function(CC,EC){
colMeans_CC <- colMeans(CC)
colMeans_EC <- colMeans(EC)
tmp <- unlist( lapply( 1:ncol(CC), function(i, colMeans_CC, colMeans_EC ){
return( ( colMeans_CC[i] - colMeans_EC[i] ) / colMeans_EC[i] )
}, colMeans_CC, colMeans_EC ) )
MD_index <- (sum(tmp) / ncol(CC)) * 100

return(MD_index)
}
