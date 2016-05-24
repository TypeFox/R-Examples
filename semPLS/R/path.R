# Used in 'read.splsm' to replace the IDs by their name.
path <-
function(variables, connections){
  ret <- matrix(NA, nrow=nrow(connections), ncol=2)
  colnames(ret) <- c("source", "target")
  for (i in 1:nrow(connections)) {
    index <- which(variables[,2]==connections[i,1])
    ret[i,1] <- as.character(variables[index,1])
    index <- which(variables[,2]==connections[i,2])
    ret[i,2] <- as.character(variables[index,1])
  }
  return(ret)
}
