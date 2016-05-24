VE.EB.SYG.Total.Hajek <- function(VecY.s, VecPk.s, MatPkl.s, N, VecAlpha.s = rep(1, times= length(VecPk.s)))
{
  if(length(N) != 1 ){stop("Value of N must be a scalar, i.e. a vector of length 1.")                                    }
  if(N%%1 != 0      ){stop("N must be an integer or a double-precision scalar with zero-valued fractional part.")        }
  DoubleN            <- as.double(N)
  VE.EB.SYG.Mean.Hajek(VecY.s, VecPk.s, MatPkl.s, VecAlpha.s) * DoubleN * DoubleN
}
