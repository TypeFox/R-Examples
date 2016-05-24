mri1 <- function(U, V, outgoing.name){
  if (missing(outgoing.name)) {UV <- cbind(U, V)} else {UV <- cbind(U, V)[V != outgoing.name, ]}
  UV <- table(UV[,1], UV[,2])
  UV. <- table(U, V)
  S <- sum(UV*(UV - 1))/2
  I <- sum(UV.*(UV. - 1))/2 + 1/2 * (sum(rowSums(UV.)**2) - sum(UV.**2))
  S/I
}
