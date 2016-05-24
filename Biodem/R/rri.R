rri <- function(x){
  somme <- colSums(x)
  x1 <- uri(x)
  dg <- diag(x1)
  conta <- matrix(0, nrow = nrow(x1), ncol = ncol(x1))
  preI1 <- rep(0, length(dg))
  for(i in 1:length(dg)){
    preI1[i] <- (dg[i]*somme[i]*(somme[i]-1))
    I1 <- sum(preI1)
  }

  for(i in 1:nrow(x1)){
    for(j in 1:ncol(x1)){
      if(i!=j)
        conta[i,j]<- (x1[i,j]*somme[i]*somme[j])
      else
        conta[i,j]=0}
    conta
  }
  preI2 <- rowSums(conta)
  I2 <- sum(preI2)

  R <- (I1 + I2)/(sum(somme) * (sum(somme)-1))
}
