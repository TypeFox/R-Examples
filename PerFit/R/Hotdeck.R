# Hotdeck imputation (dichotomous, polytomous):
HD.imputation <- function(matrix, save.matImp, ip, ability)
{
  N <- dim(matrix)[1]; I <- dim(matrix)[2]
  matrix.imp     <- matrix
  position.NA    <- is.na(matrix)
  recipients     <- which(rowSums(position.NA) > 0)
  N.recipients   <- length(recipients)
  donors         <- (1:N)[-recipients]
  N.donors       <- length(donors)
  donors.hlp.mat <- t(matrix[donors,])
  # 
  for (i in 1:N.recipients)
  {
    rcp           <- recipients[i]
    distance      <- rowSums(abs(t(matrix[rcp, ] - donors.hlp.mat)), na.rm = TRUE)
    min.distance  <- which(distance == min(distance))
    closest.donor <- if (length(min.distance) == 1) {min.distance} else {sample(min.distance, size = 1)}
    matrix.imp[rcp, position.NA[rcp, ]] <- matrix[donors[closest.donor], position.NA[rcp, ]]
  }
  # 
  if (save.matImp == TRUE)
  {
    write.matrix(matrix.imp, file="Datamatrix_imputted.txt", sep=" ")
  }
  return(list(matrix.imp, ip, ability, 1))
}