Est.EmpDistFunc.NHT <- function(VecY.s, VecPk.s, N, t)
{
  if(length(N) != 1                  ){stop("Value of N must be a scalar, i.e. a vector of length 1.")                            }
  if(N%%1 != 0                       ){stop("N must be an integer or a double-precision scalar with zero-valued fractional part.")}
  if(length(t) != 1                  ){stop("Value of t must be a scalar, i.e. a vector of length 1.")                            }
  if(t%%1 != 0                       ){stop("t must be an integer or a double-precision scalar.")                                 }
  if(! is.vector(VecY.s)             ){stop("VecY.s must be a vector.")                                                           }
  if(! is.vector(VecPk.s)            ){stop("VecPk.s must be a vector.")                                                          }
  if(any(is.na(VecPk.s))             ){stop("There are missing values in VecPk.s.")                                               }
  if(any(VecPk.s<=0|VecPk.s>1)       ){stop("There are invalid values in VecPk.s.")                                               }
  if(any(is.na(VecY.s))              ){stop("There are missing values in VecY.s.")                                                }
  n                                   <- length(VecY.s)
  if(n != as.integer(length(VecPk.s))){stop("The lengths of VecY.s and VecPk.s are different.")                                   }
  .C("Est_EmpDistFunc_NHT",
      as.double(VecY.s),
      as.double(VecPk.s),
      as.integer(n), 
      as.double(N),
      as.double(t),
      PointEst = double(1),
      PACKAGE = "samplingVarEst")$PointEst
}
