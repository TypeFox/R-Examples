similarity_hica <-
function(X,dim.subset=512){
 if (dim.subset<1){
  stop("dim.subset must be an integer positive value")
 }
 n <- dim(X)[1]
 p <- dim(X)[2]
 stat_value <- matrix(0,nrow=((p+1)-1)*((p)-1)/2,ncol=3)
 if (dim.subset>=n){
  sam=1:n
 } else {
  sam=sample(1:n,dim.subset)
 }
 count <- 1
 for (i in 1:(p-1)){
  for (j in (i+1):p){
   stat_value[count,1] <- dcor(X[sam,i],X[sam,j])
   stat_value[count,2] <- i
   stat_value[count,3] <- j
   count=count+1
  }
 }
 list(similarity_matrix=stat_value, subset=sam)
}
