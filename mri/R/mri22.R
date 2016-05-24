mri22 <- function(U, V, incoming.name){
  if (missing(incoming.name)) {UV. <- cbind(U, V)} else {UV. <- cbind(U, V)[U != incoming.name, ]}
  UV. <- table(UV.[,1], UV.[,2])
  S <- sum(UV. * (UV. - 1))/2
  I <- (sum(UV. * (UV. - 1))/2) + (sum(rowSums(UV.)**2) - sum(UV.**2))/2
  S/I
}
