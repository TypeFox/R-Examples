# Parametric model imputation
# Item scores are generated from Bernoulli distributions, with probabilities estimated by means of parametric IRT models
#     (1PLM, 2PLM, 3PLM):
PModel.imputation <- function(matrix, save.matImp, ip, model, ability, method, mu, sigma)
{
  N <- dim(matrix)[1]; I <- dim(matrix)[2]
  
  # Estimate item parameters if not provided (using 'irtoys'):
  ip <- estIP(matrix, ip, model)
  
  # Estimate ability parameters if not provided (using 'irtoys'):
  ability <- estAb(matrix, ip, ability, method, mu, sigma)
  
  A   <- ip[, 1]; B <- ip[, 2]; C <- ip[, 3]
  P   <- do.call(cbind, lapply(1:I, function (x) {C[x] + (1 - C[x]) / (1 + exp(-A[x] * (ability - B[x])))}))
  # 
  matrix.imp  <- matrix
  position.NA <- which(is.na(matrix) == 1, arr.ind = TRUE)
  P.NA        <- P[position.NA]
  matrix.imp[position.NA] <- rbinom(length(P.NA), 1, P.NA)
  # 
  if (save.matImp == TRUE)
  {
    write.matrix(matrix.imp, file="Datamatrix_imputted.txt", sep=" ")
  }
  return(list(matrix.imp, ip, ability, 1))
}
