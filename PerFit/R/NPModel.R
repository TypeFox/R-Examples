# Nonparametric model imputation (dichotomous)
# Similar to the hotdeck imputation, but item scores are generated from Bernoulli distributions, 
#    with probabilities defined by donors with similar total score than the recipient (based on all items except the NAs):
NPModel.imputation <- function(matrix, save.matImp, ip, ability)
{
  N <- dim(matrix)[1]; I <- dim(matrix)[2]
  matrix.imp   <- matrix
  position.NA  <- is.na(matrix)
  recipients   <- which(rowSums(position.NA) > 0)
  N.recipients <- length(recipients)
  donors       <- (1:N)[-recipients]
  N.donors     <- length(donors)
  # 
  vect.NC      <- rowSums(matrix, na.rm = TRUE)
  # 
  for (i in 1:N.recipients)
  {
    rcp       <- recipients[i]
    rcp.noNA  <- (1:I)[!position.NA[rcp, ]]
    rcp.NC    <- vect.NC[rcp]
    donors.NC <- rowSums(matrix[donors, rcp.noNA])
    mar       <- 0
    ctrl      <- 0
    while (ctrl == 0)
    {
      closest.donors <- (abs(donors.NC - rcp.NC) <= mar)
      if (sum(closest.donors) > 0)
      {
        ctrl <- 1
      } else
      {
        mar <- mar+1
      }
    }
    matrix.imp[rcp, position.NA[rcp, ]] <- rbinom(sum(position.NA[rcp, ]), 1, 
                                                  colMeans(matrix[donors[closest.donors], position.NA[rcp, ], drop = FALSE]))
  }
  # 
  if (save.matImp == TRUE)
  {
    write.matrix(matrix.imp, file="Datamatrix_imputted.txt", sep=" ")
  }
  return(list(matrix.imp, ip, ability, 1))
}