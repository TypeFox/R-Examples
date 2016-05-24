# MissingValues(): Dealing with missing values ----
MissingValues <- function(matrix, NAs, Save.MatImp, IP, ParModel, Ability, Method, mu, sigma)
{
  N <- dim(matrix)[1]; I <- dim(matrix)[2]
  if (sum(is.na(matrix)) > 0)
  {
    lst <- switch(
      NAs,
      Pairwise = list(matrix, IP, Ability, 1),
      Hotdeck = HD.imputation(matrix, Save.MatImp, IP, Ability),
      NPModel = NPModel.imputation(matrix, Save.MatImp, IP, Ability),
      PModel  = {
        # Sanity check - IP model:
        Sanity.IPm(ParModel)
        # Sanity check - Ability method:
        Sanity.Abm(Method)
        # 
        PModel.imputation(matrix, Save.MatImp, 
                                      IP, ParModel, Ability, Method, mu, sigma)
      }
    )
  } else
  {
    if (Save.MatImp == TRUE)
    {
      write.matrix(matrix, file="Datamatrix_original.txt", sep=" ")
    }
    lst <- list(matrix, IP, Ability, 0)
  }
  return(lst)
}

# MissingValues.poly(): Dealing with missing values (polytomous) ----
MissingValues.poly <- function(matrix, Ncat, NAs, Save.MatImp, IP, ParModel, Ability, Method)
{
  N <- dim(matrix)[1]; I <- dim(matrix)[2]
  if (sum(is.na(matrix)) > 0)
  {
    lst <- switch(
      NAs,
      Pairwise = list(matrix, IP, Ability, 1),
      Hotdeck = HD.imputation(matrix, Save.MatImp, IP, Ability),
      NPModel = NPModel.imputation.poly(matrix, Ncat, Save.MatImp, IP, Ability),
      PModel  = {
        PModel.imputation.poly(matrix, Ncat, Save.MatImp,IP, ParModel, Ability, Method)
      }
    )
  } else
  {
    if (Save.MatImp == TRUE)
    {
      write.matrix(matrix, file="Datamatrix_original.txt", sep=" ")
    }
    lst <- list(matrix, IP, Ability, 0)
  }
  lst
}
