dgambin <-
function(alpha, maxoctave)
{
  # calculates the 'fitness' distribution of species for a given alpha
  qG99 = qgamma(0.99, alpha, 1) / 100
  b1 = 0:99 * qG99
  b2 = 1:100 * qG99
  Gj <- (pgamma(b2, alpha, 1) - pgamma(b1, alpha, 1)) / 0.99
  
  # a function to give the probability distribution of the GamBin in 
  Pk <- function(k, octmax) # k is in 0:(nOct-1)
    return(sum(choose(octmax, k) * (1:100/100)^k * (1 - 1:100/100)^(octmax - k) * Gj)) 
  
  # applies Pk to each octave:
  ret <- sapply(0:maxoctave, Pk, octmax = maxoctave)
  names(ret) <- 0:maxoctave
  ret
}
