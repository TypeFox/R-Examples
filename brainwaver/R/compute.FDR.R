compute.FDR <- function(pvalue.vec,q)
{

  if (min(pvalue.vec) < 0 || max(pvalue.vec) > 1) { stop("p-values not in valid range") }

  probs <- sort(pvalue.vec) 
  l <- length(pvalue.vec)	


  correct <- sum(1/c(1:l))

  fdr <- c(1:l)/ l * (q / correct)
  sig <- array(0, c(l))

  for(i in 1:l)
  {
    if(probs[i] <= fdr[i]) sig[i] <- 1
  }

  sig <- sig * c(1:l)
  maxsig <- max(sig)
  if(maxsig==0){
    thrprob <- 0
  }else{
    thrprob <- probs[maxsig]
  }

  return(thrprob)
}