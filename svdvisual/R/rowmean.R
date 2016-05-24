rowmean <-
function(x){
   rowvec <- rowMeans(x)
   ncol <- dim(x)[2]
   ones <- rep(1,ncol)
   rowvec %*% t(ones)
}
