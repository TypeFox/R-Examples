LOGICOIL <-
function(id, seq, reg, plot.result=FALSE)
{
	# priors obtained from LOGICOIL training set
	prior <- c(0.6331,0.237,0.053,0.0769)
  
  # compute test scores
  cat("Estimating oligomeric state of coiled-coil sequences\n")
  prob.oligo <- EstimateProbability(id,
    seq,
    reg,
    pssm,
    LOGICOILfit,
    Model_Parameters)
  
  # Plot results of LOGICOIL predictions
  if(plot.result==TRUE)
  {
    plot_LOGICOIL(prob.oligo, id)
  }
  
  # compute final summary probability for whole sequence
  score <- apply(prob.oligo, 2, mean) / prior
	
  # return score fo input sequence
  cat("Analysis complete!\n")
	return(score)
	warnings()
}
