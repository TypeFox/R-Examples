compareAgreement <- function(V, n=500, e=0.01, N=500, pos=FALSE) {
  #            V = vector
  #            n = number of values to sample in error.agreement
  #            e = proportion of errors to sample in error.agreement
  #            N = number of replications for mean and standard deviation
  #          pos = positions

  # if position argument is not provided, then use the values that occur in the vector as positions
  l <- length(pos)
  if (pos[1] == FALSE & pos[l] == FALSE) pos <- as.numeric(as.data.frame(table(V))[,1])

  # input validation: enough values is pos to calculate agreement?
  if (length(pos) < 3) {
    warning("Too few positions to calculate agreement.")
    return(NA)
  }
  a <- agreement(collapse(V, pos=pos))
  z <- replicate(N, agreementError(V, n=n, e=e, pos=pos))
  return(list(agreement=a, mean=mean(z), sd=sd(z)))
}
