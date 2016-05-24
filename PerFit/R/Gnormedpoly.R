########################################################################################
########################################################################################
# Polytomous items: Number of Guttman errors
# Gpoly reduces to G for 0-1 data (Ncat=2) (currently Gpoly needs G to be loaded)
# Gnormed.poly reduces to Gnormed for 0-1 data (Ncat=2) (currently Gnormed.poly needs G to be loaded)
########################################################################################
########################################################################################

Gnormed.poly <- function(matrix, Ncat,
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
  # Numerator:
  probs.ISD      <- matrix(NA, nrow = I, ncol = M)
  for (m in 1:M) {probs.ISD[, m] <- colMeans(matrix >= m, na.rm = TRUE)}
  f.scoresISD    <- function (x) {if (is.na(x)) {rep(NA, M)} else {c(rep(1, x), rep(0, M - x))}}
  matrix.ISD     <- matrix(unlist(lapply(t(matrix), f.scoresISD)), byrow = TRUE, nrow = N)
  probs.ISD.vect <- as.vector(t(probs.ISD))
  matrix.ISD.ord <- matrix.ISD[, order(probs.ISD.vect, decreasing = TRUE)]
  num            <- G(matrix.ISD.ord)$PFscores[, 1]
  # Denominator: 
  NC        <- rowSums(matrix, na.rm = TRUE)
  res <- rep(NA, N)
  # A vectorial solution is possible without NAs (until PerFit 1.3.1).
  # With NAs a respondent-by-respondent approach is needed.
  # This is here implemented in two steps: 
  #    (I)  Respondents without NAs (in one go).
  #    (II) Respondents with NAs (one by one).
  # Step (I):
  rows.I <- rowSums(is.na(matrix)) == 0
  if (length(rows.I) > 0)
  {
    ranks.ISD <- matrix(rank(I * M - probs.ISD.vect, ties.method = "first"), nrow = I, byrow = TRUE)
    if (Ncat>2) {cumranks.ISD <- cbind(rep(0, I), t(apply(ranks.ISD, 1, cumsum)))}
    if (Ncat==2) {cumranks.ISD <- cbind(rep(0, I), ranks.ISD)}
    V        <- matrix(rep(cumranks.ISD[1, ], Ncat), ncol = Ncat, byrow = FALSE)
    add.term <- matrix(rep(cumranks.ISD[2, ], nrow(V)), ncol = Ncat, byrow = TRUE)
    V        <- V + add.term
    T        <- sapply(2:sum(dim(V)), function(x) {max(V[col(V)+row(V) == x])})
    for (i in 3:I) 
    {
      V        <- matrix(rep(T, Ncat), ncol=Ncat, byrow=FALSE)
      add.term <- matrix(rep(cumranks.ISD[i,], nrow(V)), ncol=Ncat, byrow=TRUE)
      V        <- V + add.term
      T        <- sapply(2:sum(dim(V)), function(x) {max(V[col(V)+row(V) == x])})
    }
    maxG       <- T - sapply(0:(I * M), function(x) {.5 * x * (x+1)})
    # 
    res[rows.I] <- num[rows.I] / maxG[NC[rows.I] + 1]
  }
  # Step (II):
  rows.II <- rowSums(is.na(matrix)) > 0
  if (length(rows.II) > 0)
  {
    den.II <- apply(matrix[rows.II, ], 1, function(vec)
    {
      if (sum(is.na(vec)) == I)
      {
        return(0)
      } else
      {
        Items               <- (1:I)[!is.na(vec)]
        I.noNA              <- length(Items)
        probs.ISD.vect.noNA <- as.vector(t(probs.ISD[Items, ]))
        ranks.ISD           <- matrix(rank(I.noNA * M - probs.ISD.vect.noNA, ties.method = "first"), nrow = I.noNA, byrow = TRUE)
        if (Ncat>2) {cumranks.ISD <- cbind(rep(0, I.noNA), t(apply(ranks.ISD, 1, cumsum)))}
        if (Ncat==2) {cumranks.ISD <- cbind(rep(0, I.noNA), ranks.ISD)}
        V        <- matrix(rep(cumranks.ISD[1, ], Ncat), ncol = Ncat, byrow = FALSE)
        add.term <- matrix(rep(cumranks.ISD[2, ], nrow(V)), ncol = Ncat, byrow = TRUE)
        V        <- V + add.term
        T        <- sapply(2:sum(dim(V)), function(x) {max(V[col(V)+row(V) == x])})
        for (i in 3:I.noNA) 
        {
          V        <- matrix(rep(T, Ncat), ncol = Ncat, byrow = FALSE)
          add.term <- matrix(rep(cumranks.ISD[i,], nrow(V)), ncol = Ncat, byrow = TRUE)
          V        <- V + add.term
          T        <- sapply(2:sum(dim(V)), function(x) {max(V[col(V)+row(V) == x])})
        }
        maxG       <- T - sapply(0:(I.noNA * M), function(x) {.5 * x * (x+1)})
        return(maxG[sum(vec, na.rm = TRUE) + 1])
      }
    })
    res[rows.II] <- num[rows.II] / den.II
  }
  # 
  res[is.nan(res)] <- 0 
  # Export results:
  export.res.NP(matrix, N, res, "Gnormed.poly", vector("list", 5) , Ncat=Ncat, NA.method, 
                IRT.PModel, res.NA[[2]], Ability.PModel, res.NA[[3]], IP.NA, Ability.NA, res.NA[[4]])
  
}