########################################################################################
########################################################################################
# lzp
########################################################################################
########################################################################################
lzpoly <- function(matrix, Ncat,
                   NA.method="Pairwise", Save.MatImp=FALSE, 
                   IP=NULL, IRT.PModel="GRM", Ability=NULL, Ability.PModel="EAP")
{
  #
  matrix <- as.matrix(matrix)
  N      <- dim(matrix)[1]; I <- dim(matrix)[2]; M <- Ncat-1
  IP.NA  <- is.null(IP); Ability.NA  <- is.null(Ability)
  # Sanity check - Data matrix adequacy:
  Sanity.dma.poly(matrix, N, I, M)
  # Estimate item parameters if not provided (polytomous):
  IP.res <- estIP.poly(matrix, Ncat, IP, IRT.PModel)
  IP     <- IP.res[[1]]
  IP.ltm <- IP.res[[2]] 
  # Estimate ability parameters if not provided (using 'ltm'):
  Ability <- estAb.poly(matrix, IP.ltm, Ability, Ability.PModel)
  # Dealing with missing values:
  res.NA <- MissingValues.poly(matrix, Ncat, NA.method, Save.MatImp, IP, IRT.PModel, Ability, Ability.PModel)
  matrix <- res.NA[[1]]
  # Perfect response vectors allowed.
  # Compute PFS:
  P.CRF   <- estP.CRF(I, Ncat, IRT.PModel, IP, Ability)
  # 
  idty      <- diag(Ncat)
  f.scores  <- function (x) {idty[x+1, ]}
  matrix.01 <- matrix(unlist(lapply(t(matrix), f.scores)), byrow = TRUE, nrow = N)
  # If there are answer options not chosen by any respondent then some entries in 'P.CRF' might be 0.
  # Below all corresponding logs are set from Inf to 0.
  # (Reason: They carry no information regarding aberrant response behavior).
  log.P.CRF <- log(P.CRF)
  log.P.CRF[is.infinite(log.P.CRF)] <- 0           
  #
  l0p  <- rowSums(matrix.01 * log.P.CRF, na.rm = TRUE)
  El0p <- rowSums(P.CRF * log.P.CRF)
  # Variance (two equivalent options up to time efficiency):
  if (I*Ncat < 300)
  {
    ones.block <- vector("list", I)
    ones.block[1:I] <- list(matrix(rep(1, Ncat^2), nrow = Ncat))
    to.sum <- as.matrix(do.call(bdiag, ones.block))
    V.row  <- function(vect) {
      log.vect <- log.P.CRF[vect[1],]
      vect2 <- vect[2:(dim(P.CRF)[2]+1)]
      tmp.part1 <- (vect2 %*% t(vect2)) * to.sum
      tmp.part2 <- matrix(rep(log.vect, dim(P.CRF)[2]), nrow=dim(P.CRF)[2]) * to.sum
      # 
      sum(tmp.part1 * tmp.part2 * (tmp.part2 - t(tmp.part2)))
    }
    Vl0p <- apply(cbind(1:N, P.CRF), 1, V.row)
  } else
  {
    V.row <- function(vect) {
      tot <- 0;
      for (i in 1:I) {
        log.vect  <- log.P.CRF[vect[1], ((i-1)*Ncat+1):(Ncat*i)]
        vect2     <- vect[((i-1)*Ncat+1):(Ncat*i)+1]
        tmp.part1 <- vect2 %*% t(vect2)
        tmp.part2 <- matrix(rep(log.vect, Ncat), nrow=Ncat)
        # 
        tot       <- tot + sum(tmp.part1 * tmp.part2 * (tmp.part2 - t(tmp.part2)));
      }
      tot
    }
    Vl0p <- apply(cbind(1:N,P.CRF),1,V.row)
  }
  res <- (l0p - El0p) / sqrt(Vl0p)
  res[rowSums(is.na(matrix)) == I] <- NA
  # Export results:
  export.res.P(matrix, N, res, "lzpoly", vector("list", 5) , Ncat=Ncat, NA.method, 
               IRT.PModel, res.NA[[2]], Ability.PModel, res.NA[[3]], IP.NA, Ability.NA, res.NA[[4]])
}