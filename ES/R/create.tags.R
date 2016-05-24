create.tags <-
function(mat){
  n <- NROW(mat)
  k <- NCOL(mat)
  cnms <- colnames(mat)

  resind <- matrix(0, nrow=k*(k-1)/2, ncol=2)
  
  ind1 <- 1:n
  j0 <- 1
  for(i in 1:(k-1)){
    ind2 <- ind1 + n
    for(j in (i+1):k){
           j0 <- j0+1

      resind[j0/2, 1] <- i
      resind[j0/2, 2] <- j
      
      j0 <- j0 + 1
      ind2 <- ind2 + n
    }
    ind1 <- ind1 + n
  }
   col.means <- colMeans(mat) 
     col.norms <- sqrt(colSums(mat*mat))
  list(tags=resind,col.means=col.means, col.norms=col.norms)
}