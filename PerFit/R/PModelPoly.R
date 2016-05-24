# Parametric model imputation (polytomous)
# Item scores are generated from multinomial distributions, with probabilities estimated by means of parametric IRT models
#     (PCM, GPCM, GRM):
PModel.imputation.poly <- function(matrix, Ncat, save.matImp, ip, model, ability, method)
{
  N <- dim(matrix)[1]; I <- dim(matrix)[2]
  # 
  matrix2 <- data.frame(apply(matrix, 2, as.factor)) # eliminates item levels with no answers
  
  # Estimate item parameters if not provided (polytomous):
  ip.res <- estIP.poly(matrix, Ncat, ip, model)
  ip     <- ip.res[[1]]
  ip.ltm <- ip.res[[2]] 
  
  # Estimate ability parameters if not provided (using 'ltm'):
  ability <- estAb.poly(matrix, ip.ltm, ability, method)
  
  # Compute P.CRF:
  P.CRF <- estP.CRF(I, Ncat, model, ip, ability)
  
  # Fill in NAs:
  matrix.imp      <- matrix
  position.NA.mat <- is.na(matrix)
  resp.NA         <- which(rowSums(position.NA.mat) > 0)
  for (i in 1:length(resp.NA))
  {
    resp <- resp.NA[i]
    position.NA <- (1:I)[position.NA.mat[resp, ]]
    P.CRF.ind   <- sapply(position.NA,function(x){((x - 1) * Ncat + 1):(x * Ncat)})
    P.CRF.NA    <- matrix(P.CRF[resp, P.CRF.ind], ncol = length(position.NA))
    matrix.imp[resp, position.NA] <- 
      which(apply(P.CRF.NA, 2, function(vect) {rmultinom(1, 1, vect)}) == 1, arr.ind = TRUE)[, 1] - 1
  }
  # 
  if (save.matImp == TRUE)
  {
    write.matrix(matrix.imp, file="Datamatrix_imputted.txt", sep=" ")
  }
  return(list(matrix.imp, ip, ability, 1))
}
