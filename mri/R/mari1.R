mari1 <- function(U, V, outgoing.name, k = 500){
  R <- NULL
  for (i in 1:k) R[i] <- mri1(sample(U), sample(V), outgoing.name=outgoing.name)
  O <- mri1(U, V, outgoing.name=outgoing.name)
  E <- mean(R)
  (O - E)/(1 - E)
}
