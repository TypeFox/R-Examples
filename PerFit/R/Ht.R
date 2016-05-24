########################################################################################
########################################################################################
# Ht (Sijtsma, 1986)
########################################################################################
########################################################################################
Ht <- function(matrix, 
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
  # Sanity check - Perfect response vectors:
  part.res  <- Sanity.prv(matrix, N, I)
  NC        <- part.res$NC
  all.0s    <- part.res$all.0s
  all.1s    <- part.res$all.1s
  matrix.sv <- matrix
  matrix    <- part.res$matrix.red
  # Compute PFS: 
  singlePs  <- rowMeans(matrix, na.rm = TRUE)
  tot.score <- colSums(matrix.sv, na.rm = TRUE)
  N.red     <- dim(matrix)[1]
  I.NAs     <- I - rowSums(is.na(matrix))
  num       <- apply(matrix, 1, function(vect){cov(vect, tot.score - vect, use = "pairwise.complete.obs")}) * (I.NAs - 1) / I.NAs
  df        <- data.frame(1:N.red, singlePs, num)
  df.ord    <- df[order(df[, 2]), ]
  singlePs.ord <- df.ord[[2]]
  pos       <- which(diff(singlePs.ord, lag = 1) > 0) 
  less      <- sapply(c(pos, N.red), function(x)
  {
    sum(singlePs.ord[1:(x-1)])
  })
  if (pos[1]==1){less[1] <- 0}
  less      <- rep(less, c(pos[1], diff(c(pos, N.red), lag = 1))) * (1 - singlePs.ord)
  more      <- sapply(pos, function(x){sum(1 - singlePs.ord[(x+1):N.red])}); more <- c(more, 0)
  more      <- rep(more, c(pos[1], diff(c(pos, N.red), lag = 1))) * singlePs.ord
  den       <- less + more
  df.ord    <- data.frame(df.ord, den)
  df        <- df.ord[order(df.ord[, 1]),]
  res.red   <- df$num / df$den
  # Compute final PFS vector:
  res <- final.PFS(res.red, all.0s, all.1s, N)
  # Export results:
  export.res.NP(matrix.sv, N, res, "Ht", part.res, Ncat=2, NA.method, 
                IRT.PModel, res.NA[[2]], Ability.PModel, res.NA[[3]], IP.NA, Ability.NA, res.NA[[4]])
}
