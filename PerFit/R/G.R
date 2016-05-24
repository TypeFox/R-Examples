########################################################################################
########################################################################################
# G (van der Flier, 1977; Meijer, 1994):
########################################################################################
########################################################################################
G <- function(matrix, 
              NA.method="Pairwise", Save.MatImp=FALSE, 
              IP=NULL, IRT.PModel="2PL", Ability=NULL, Ability.PModel="ML", mu=0, sigma=1)
{
  matrix      <- as.matrix(matrix)
  N           <- dim(matrix)[1]; I <- dim(matrix)[2]
  IP.NA       <- is.null(IP); Ability.NA  <- is.null(Ability)
  # Sanity check - Data matrix adequacy:
  Sanity.dma(matrix, N, I)
  # Dealing with missing values:
  res.NA <- MissingValues(matrix, NA.method, Save.MatImp, IP, IRT.PModel, Ability, Ability.PModel, mu, sigma)
  matrix <- res.NA[[1]]
  # Perfect response vectors allowed (albeit uninformative).
  # Compute PFS:
  NC         <- rowSums(matrix, na.rm = TRUE)
  uniqueNC   <- sort(unique(NC))
  pi         <- colMeans(matrix, na.rm = TRUE)
  matrix.ord <- matrix[, order(pi, decreasing = TRUE)]
  per.row    <- function(vect)
  {
    ind.0     <- which(vect == 0)
    ind.1     <- which(vect == 1)
    all.cases <- expand.grid(ind.0, ind.1)
    sum((all.cases[, 2] - all.cases[, 1]) > 0)
  }
  res <- apply(matrix.ord, 1, per.row)
  res[rowSums(is.na(matrix)) == I] <- NA
  # Export results:
  export.res.NP(matrix, N, res, "G", vector("list", 5), Ncat=2, NA.method, 
                IRT.PModel, res.NA[[2]], Ability.PModel, res.NA[[3]], IP.NA, Ability.NA, res.NA[[4]])
}
