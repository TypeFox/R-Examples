getTask <-
function(n){
  tasks <- matrix(0, nrow = (n*n + n)/2, ncol = 2)
  iN <- seq_len(n)
  tasks[, 1] <- rep(iN, sort(iN, decreasing = T))
  tasks[, 2] <- unlist(lapply(iN, function(x) iN[x:n]))
  colnames(tasks) <- paste("Var", 1:2)
  return(tasks)
}
