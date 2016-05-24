scorer2 <-
function(id, seq, reg, pssm, delta=1)
{
  # load training PSSM
  #cat("Loading training PSSM\n")
  
  # Define levels of variables in model ###
  var <- list(
      amino = c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y","X"),
      register = letters[1:7])
 
  # set delta parameter
  #delta <- 1

  # compute test scores
  cat("Estimating oligomeric state of coiled-coil sequences\n")
  score <- EstimateProbability(id, seq, reg, pssm, var, delta)

  cat("Analysis complete!\n")
  out <- data.frame(ID = id, score = score)
  return(out)
}
