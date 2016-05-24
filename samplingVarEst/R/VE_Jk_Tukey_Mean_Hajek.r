VE.Jk.Tukey.Mean.Hajek <- function(VecY.s, VecPk.s, N, FPC= TRUE)
{
  VE.Jk.Tukey.Ratio(VecY.s, as.double(rep(1, times= length(VecY.s))), VecPk.s, N, FPC)
}
