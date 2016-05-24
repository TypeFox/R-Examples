nnWCM <- function(matrix, num_cat) {
  
  #Function computes the non-normalized within-cluster mutability (based on the Gini coefficient)
  dim <- dim(matrix)
  #number of variables
  m <- dim[3]
  #number of clusters
  k <- dim[2]
  #number of elements
  n <- sum(matrix[,,1])
  
  var <- vector(mode="numeric", length=k)
  
  for (g in 1:k) {
    gini <- vector(mode="numeric", length=m)
    for (l in 1:m) {
      step <- 0
      step <- matrix[,g,l]
      
      K <- num_cat[l]
      cluster <- sum(step)
      
      temp <- 0
      for (u in 1:K) {
        temp[u] <- (step[u]/cluster)^2
      }
      gini[l] <- (1-sum(temp))
    }
    var[g] <- sum(gini)*cluster/n
  }
  sum(var)
}