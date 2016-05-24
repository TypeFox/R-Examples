POS <- function(ES, Core, Y){
  Y <- as.factor(Y)
  levels(Y) = c(1,2)
  # total core length
  t.Length = pmax(Core[,2], Core[,4]) - pmin(Core[,1], Core[,3])
  # overlap length
  b <- pmin(Core[,2], Core[,4])
  a <- pmax(Core[,1], Core[,3])
  o.Length = b - a
  o.Length[o.Length < 0] <- 0
  # Number of inlier (non-outlier) observations
  L <- ES
  L[, Y == 1] = (ES[, Y == 1] >= Core[,1] & ES[, Y == 1] <= Core[,2])
  L[, Y == 2] = (ES[, Y == 2] >= Core[,3] & ES[, Y == 2] <= Core[,4])
  n.inlier = apply(L, 1, sum)
  # Number of overlapping observations
  L[o.Length > 0, ] = L[o.Length > 0, ] & (ES[o.Length > 0,] <= b[o.Length > 0]) & (ES[o.Length > 0,] >= a[o.Length > 0])
  L[o.Length <= 0, ] = FALSE
  n.ovrlap1 <- apply(L[, Y == 1], 1, sum)
  n.ovrlap2 <- apply(L[, Y == 2], 1, sum)
  n.overlap = n.ovrlap1 + n.ovrlap2
  # POS
  POS = 4*(o.Length/t.Length)*(n.ovrlap1*n.ovrlap2/n.inlier)
  POS[POS != 0] = POS[POS != 0]/n.overlap[POS != 0]   
  POS[is.na(POS)] = 1   #in Some Proteomic data, all observations may have a specific value for some genes ==> So, t.Length=0 ==> POS=NA
  return(POS)
}