VE.EB.SYG.Mean.Hajek <- function(VecY.s, VecPk.s, MatPkl.s, VecAlpha.s = rep(1, times= length(VecPk.s)))
{
  VE.EB.SYG.Ratio(VecY.s, as.double(rep(1, times= length(VecY.s))), VecPk.s, MatPkl.s, VecAlpha.s)
}
