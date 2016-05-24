uri <- function(x){
  somme <- colSums(x)
  M <- matrix(data=0,nrow=length(somme),ncol=length(somme))
  for(i in 1:length(somme)){
    for(j in 1:length(somme)){
      if(i==j)
        M[i,j] <- (sum(x[,i]*(x[,i]-1)))/(somme[i]*(somme[i]-1))
      else
        M[i,j] <- (sum(x[,i]*x[,j]))/(somme[i]*somme[j])
    }
  }
  M
}
