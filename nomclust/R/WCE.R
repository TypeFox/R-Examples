WCE <- function(matrix, num_cat) {
  #Function computes the normalized within-cluster entropy
  dim <- dim(matrix)
  #number of variables
  m <- dim[3]
  #number of clusters
  k <- dim[2]
  #number of elements
  n <- sum(matrix[,,1])
  
  var <- vector(mode="numeric", length=k)
  for (g in 1:k) {
    entropy_norm <- vector(mode="numeric", length=m)
    for (l in 1:m) {
      step <- 0
      step <- matrix[,g,l]
      
      K <- num_cat[l]
      cluster <- sum(step)
      
      temp <- 0
      for (u in 1:K) {
        if (step[u] == 0) {
          temp[u] <- 0
        }
        else {
          temp[u] <- ((step[u]/cluster)%*%log(step[u]/cluster))
        }
      }
      entropy_norm[l] <- -sum(temp)/log(K)
    }
    var[g] <- sum(entropy_norm)*cluster/(m*n)
  }
  sum(var)
}