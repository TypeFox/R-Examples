VE.Jk.EB.SW2.Total.Hajek <- function(VecY.s, VecPk.s, nII, VecPi.s, VecCluLab.s, VecCluSize.s, N)
{
  if(length(N) != 1 ){stop("Value of N must be a scalar, i.e. a vector of length 1.")                                    }
  if(N%%1 != 0      ){stop("N must be an integer or a double-precision scalar with zero-valued fractional part.")        }
  DoubleN            <- as.double(N)
  VE.Jk.EB.SW2.Mean.Hajek(VecY.s, VecPk.s, nII, VecPi.s, VecCluLab.s, VecCluSize.s) * DoubleN * DoubleN
}
