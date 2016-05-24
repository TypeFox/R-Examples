HMMcopynumber <-
function(obs, transitionmatrix, emissionmatrix, includeZeroState, windowlength, chrlength)
{
  states <- rownames(emissionmatrix)
  cat("Attempting to create Viterbi matrix\n")
  v <- makeViterbimat(obs, transitionmatrix, emissionmatrix, includeZeroState)
  
  # Go through each of the rows of the matrix v (where each row represents
  # a window in the DNA sequence), and find out which column has the
  # maximum value for that row.
  
  mostprobablestatepath <- apply(v, 1, function(x) which.max(x))
  
  # Create table of breakpoints and most likely states
  
  results <- data.frame(Startpos=numeric(0), Endpos=numeric(0), State=numeric(0), stringsAsFactors=F)
  
  prevobs <- obs[1]
  prevmostprobablestate <- mostprobablestatepath[1]
  prevmostprobablestatename <- states[prevmostprobablestate]
  startpos <- 1
  for (i in 2:length(obs))
  {
    observation <- obs[i]
    mostprobablestate <- mostprobablestatepath[i]
    mostprobablestatename <- states[mostprobablestate]
    if (mostprobablestatename != prevmostprobablestatename)
    {
      cat("Positions",startpos,"-",(i-1)*windowlength, ": Most probable state = ", prevmostprobablestatename,"\n")
      results[nrow(results)+1,] <- c(startpos, (i-1)*windowlength, prevmostprobablestatename)
      startpos <- (i-1)*windowlength+1
    }
    prevobs <- observation
    prevmostprobablestatename <- mostprobablestatename
  }
  cat("Positions",startpos,"-",chrlength, ": Most probable state = ", prevmostprobablestatename,"\n")
  results[nrow(results)+1,] <- c(startpos, chrlength, prevmostprobablestatename)
  
  return(results)
}
