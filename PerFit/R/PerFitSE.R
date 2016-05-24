# Jackknife procedure to estimate SE for PFSs:
PerFit.SE <- function(x)
{
  dico.PFS <- c("Cstar", "C.Sato", "U3", "ZU3", "G", "Gnormed", "D.KB", "r.pbis", "NCI", "Ht", "A.KB", "E.KB", "lz", "lzstar")
  poly.PFS <- c("Gpoly", "Gnormed.poly", "U3poly", "lzpoly")
  # 
  matrix <- x$Matrix      # 'matrix' is NA-free
  N      <- dim(matrix)[1]; I <- dim(matrix)[2] 
  PFS    <- x$PFStatistic
  JK.mat <- matrix(NA, nrow = N, ncol = I)
  if (PFS == "lzpoly")
  {
    for (it in 1:I)
    {
      JK.mat[, it] <- do.call(PFS, c( list(matrix[, -it]), list(Ncat = x$Ncat, IP = NULL, Ability = NULL) ))$PFscores[, 1]
    }
  } else
  {
    if (PFS %in% dico.PFS)
    {
      for (it in 1:I)
      {
        JK.mat[, it] <- do.call(PFS, c( list(matrix[, -it]), list(IP = x$IP[-it, ], Ability = NULL) ))$PFscores[, 1]
      }
    }
    if (PFS %in% poly.PFS)
    {
      for (it in 1:I)
      {
        JK.mat[, it] <- do.call(PFS, c( list(matrix[, -it]), list(Ncat = x$Ncat, IP=x$IP[-it, ], Ability = NULL) ))$PFscores[, 1]
      }
    }
  }
  I.NAs  <- I - rowSums(is.na(JK.mat))
  SE     <- sqrt( ((I.NAs - 1)/I.NAs ) * rowSums((JK.mat - rowMeans(JK.mat, na.rm = TRUE))^2, na.rm = TRUE) )
  PFS.NA <- c(x$ID.all0s, x$ID.all1s)
  if (!is.null(PFS.NA))
  {
    SE[PFS.NA] <- NA
  }
  return(cbind(PFscores = x$PFscores[, 1], PFscores.SE = round(SE, 4)))
}
