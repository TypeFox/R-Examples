########################################################################################
########################################################################################
# lzstar
########################################################################################
########################################################################################
lzstar <- function(matrix, 
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
  P   <- do.call(cbind, lapply(1:I, function (x) {C[x] + (1 - C[x]) / (1 + exp(-A[x] * (Ability - B[x])))}))
  Q   <- 1-P
  d1P <- do.call(cbind, lapply(1:I, function (x){
    (1 - C[x]) * A[x] * exp(A[x] * (Ability - B[x])) / (1 + exp(A[x] * (Ability - B[x])))^2}))
  d2P <- do.call(cbind, lapply(1:I, function (x){
    (1 - C[x]) * (A[x]^2) * exp(A[x] * (Ability - B[x])) * (1 - exp(A[x] * (Ability - B[x]))) / (1 + exp(A[x] * (Ability - B[x])))^3}))
  ri  <- d1P / (P * Q)
  r0  <- switch(Ability.PModel,
               ML = 0,
               BM = (mu - Ability) / (sigma^2),
               WL = rowSums((d1P * d2P) / (P * Q)) / (2 * rowSums((d1P^2) / (P * Q))))
  wi       <- log(P/Q)
  Wn       <- rowSums((matrix - P)*wi, na.rm = TRUE)
  sigma2n  <- rowSums((wi^2) * P * Q) / I
  cn       <- rowSums(d1P * wi) / rowSums(d1P * ri)
  wi.tilde <- wi - matrix(rep(cn, I), nrow = N) * ri
  tau2n    <- rowSums((wi.tilde^2) * P * Q) / I
  EWn      <- -cn * r0
  VWn      <- I * tau2n
  res      <- as.vector(round((Wn - EWn) / sqrt(VWn), 4))
  # Export results:
  export.res.P(matrix, N, res, "lzstar", vector("list", 5) , Ncat=2, NA.method, 
               IRT.PModel, res.NA[[2]], Ability.PModel, res.NA[[3]], IP.NA, Ability.NA, res.NA[[4]])
}