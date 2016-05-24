mari22 <- function(U, V, incoming.name, k = 500){
  R <- NULL
  for (i in 1:k) R[i] <- mri22(sample(U), sample(V), incoming.name = incoming.name)
  O <- mri22(U, V, incoming.name = incoming.name)
  E <- mean(R)
  (O - E)/(1 - E)
}
