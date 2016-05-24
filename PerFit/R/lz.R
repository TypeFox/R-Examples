########################################################################################
########################################################################################
# lz
########################################################################################
########################################################################################
lz <- function(matrix, 
                NA.method="Pairwise", Save.MatImp=FALSE, 
                IP=NULL, IRT.PModel="2PL", Ability=NULL, Ability.PModel="ML", mu=0, sigma=1)
{
  matrix      <- as.matrix(matrix)
  N           <- dim(matrix)[1]; I <- dim(matrix)[2]
  IP.NA       <- is.null(IP); Ability.NA  <- is.null(Ability)
  # Sanity check - Data matrix adequacy:
  Sanity.dma(matrix, N, I)
  # Estimate item parameters if not provided (using 'irtoys'):
  IP <- estIP(matrix, IP, IRT.PModel)
  # Estimate ability parameters if not provided (using 'irtoys'):
  Ability <- estAb(matrix, IP, Ability, Ability.PModel, mu, sigma)
  # Dealing with missing values:
  res.NA <- MissingValues(matrix, NA.method, Save.MatImp, IP, IRT.PModel, Ability, Ability.PModel, mu, sigma)
  matrix <- res.NA[[1]]
  # Perfect response vectors allowed.
  # Compute PFS: 
  A   <- IP[, 1]; B <- IP[, 2]; C <- IP[, 3]
  P   <- do.call(cbind, lapply(1:I, function (x) {C[x] + (1 - C[x]) / (1 + exp(-A[x] * (Ability - B[x])))})); Q <- 1-P
  l0  <- rowSums(matrix * log(P) + (1 - matrix) * log(Q), na.rm = TRUE)
  El0 <- rowSums(P * log(P) + Q * log(Q), na.rm = TRUE)
  Vl0 <- rowSums(P * Q * (log(P / Q))^2, na.rm = TRUE)
  res <- as.vector(round((l0 - El0) / sqrt(Vl0),4))
  # Rows of 'matrix' full of NAs -> lz = NA
  res[is.nan(res)] <- NA
  # Export results:
  export.res.P(matrix, N, res, "lz", vector("list", 5) , Ncat=2, NA.method, 
                IRT.PModel, res.NA[[2]], Ability.PModel, res.NA[[3]], IP.NA, Ability.NA, res.NA[[4]])
}