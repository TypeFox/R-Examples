RetrainScorer2 <-
function(seq, reg, type)
{
  training.data <- cbind(seq, reg, type)
  colnames(training.data) <- c("sequence", "register", "type")
  # Define levels of variables in model ###
  var <- list(
      amino = c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y","X"),
      register = letters[1:7])
  
  training.count <- nrow(training.data)
  
  # compute profile scoring matrices
  cat("Computing PSSMs\n")
  pssm <- CreatePssm(training.data, var)

  return(pssm)
}
