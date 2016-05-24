lasker <- function(x){
  somme <- colSums(x)
  M <- matrix(data=0,nrow=length(somme),ncol=length(somme))
  rownames(M) <- dimnames(x)[[2]]
  colnames(M) <- dimnames(x)[[2]]
  for(i in 1:length(somme)){
    for(j in 1:length(somme)){
      M[i,j] <- (sum(x[,i]*x[,j]))/(2*(somme[i]*somme[j]))
    }
  }
  M
}

