########################################################################################
########################################################################################
# Polytomous items: Number of Guttman errors
# Gpoly reduces to G for 0-1 data (Ncat=2) (currently Gpoly needs G to be loaded)
# Gnormed.poly reduces to Gnormed for 0-1 data (Ncat=2) (currently Gnormed.poly needs G to be loaded)
########################################################################################
########################################################################################

Gpoly <- function(matrix, Ncat,
                  NA.method="Pairwise", Save.MatImp=FALSE, 
                  IP=NULL, IRT.PModel="GRM", Ability=NULL, Ability.PModel="EAP")
{
  matrix      <- as.matrix(matrix)
  N           <- dim(matrix)[1]; I <- dim(matrix)[2]; M <- Ncat-1
  IP.NA       <- is.null(IP); Ability.NA  <- is.null(Ability)
  # Sanity check - Data matrix adequacy:
  Sanity.dma.poly(matrix, N, I, M)
  # Dealing with missing values:
  res.NA <- MissingValues.poly(matrix, Ncat, NA.method, Save.MatImp, IP, IRT.PModel, Ability, Ability.PModel)
  matrix <- res.NA[[1]]
  # Perfect response vectors allowed (albeit uninformative).
  # Compute PFS:
  probs.ISD       <- matrix(NA, nrow = I, ncol = M)
  for (m in 1:M) {probs.ISD[, m] <- colMeans(matrix >= m, na.rm = TRUE)}
  f.scoresISD    <- function (x) {if (is.na(x)) {rep(NA, M)} else {c(rep(1, x), rep(0, M - x))}}
  matrix.ISD      <- matrix(unlist(lapply(t(matrix), f.scoresISD)), byrow = TRUE, nrow = N)
  probs.ISD.vect  <- as.vector(t(probs.ISD))
  matrix.ISD.ord  <- matrix.ISD[, order(probs.ISD.vect, decreasing = TRUE)]
  res             <- G(matrix.ISD.ord)$PFscores[, 1]
  # Export results:
  export.res.NP(matrix, N, res, "Gpoly", vector("list", 5) , Ncat=Ncat, NA.method, 
               IRT.PModel, res.NA[[2]], Ability.PModel, res.NA[[3]], IP.NA, Ability.NA, res.NA[[4]])
}
