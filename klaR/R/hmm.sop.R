hmm.sop <- function(sv, trans.matrix, prob.matrix)
{
# sv: Start value
# trans.matrix: probability matrix of transitions
# prob.matrix: p(c|x)
  rekurs <- prob.matrix
  rekurs[1,] <- prob.matrix[1, ] * trans.matrix[sv, ]
  for(z in 2:nrow(prob.matrix))
    for(j in 1:ncol(prob.matrix))
        rekurs[z,j] <- 
          prob.matrix[z, j] * trans.matrix[, j] %*% rekurs[(z-1), ]
  rekurs <- rekurs / rowSums(rekurs)
  return(rekurs)
}
