VE.Jk.Tukey.RegCo.Hajek <- function(VecY.s, VecX.s, VecPk.s, N, FPC= TRUE)
{
  if(! is.vector(VecY.s)              ){stop("VecY.s must be a vector.")                                                                   }
  if(! is.vector(VecX.s)              ){stop("VecX.s must be a vector.")                                                                   }
  if(! is.vector(VecPk.s)             ){stop("VecPk.s must be a vector.")                                                                  }
  if(any(is.na(VecPk.s))              ){stop("There are missing values in VecPk.s.")                                                       }
  if(any(VecPk.s<=0|VecPk.s>1)        ){stop("There are invalid values in VecPk.s.")                                                       }
  if(any(is.na(VecY.s))               ){stop("There are missing values in VecY.s.")                                                        }
  if(any(is.na(VecX.s))               ){stop("There are missing values in VecX.s.")                                                        }
  if((length(N) != 1) & FPC           ){stop("Value of N must be a scalar, i.e. a vector of length 1.")                                    }
  if((N%%1      != 0) & FPC           ){stop("N must be an integer or a double-precision scalar with zero-valued fractional part.")        }
  n                                    <- length(VecY.s)
  if(n != length(VecPk.s)             ){stop("The lengths of VecY.s and VecPk.s are different.")                                           }
  if(n != length(VecX.s)              ){stop("The lengths of VecY.s and VecX.s are different.")                                            }
  if(length(VecPk.s) != length(VecX.s)){stop("The lengths of VecPk.s and VecX.s are different.")                                           }
  DoubleN                              <- as.double(N)
  Doublen                              <- as.double(n)
  VecEstTheta_k                        <- .C("Est_RegCo_Hajek_Excluding_All_Elements",
                                              as.double(VecY.s),
                                              as.double(VecX.s),
                                              as.double(VecPk.s),
                                              as.integer(n),
                                              VectVarEst = double(n),
                                              PACKAGE = "samplingVarEst")$VectVarEst
  EstTheta                             <- Est.RegCo.Hajek(VecY.s, VecX.s, VecPk.s)
  OUTPUT                               <- (Doublen-1)/Doublen * sum( (EstTheta - VecEstTheta_k)^2 )
  if(FPC                              ){OUTPUT <- (1-Doublen/DoubleN) * OUTPUT                                                             }
  if(OUTPUT<0                         ){warning("The variance estimate contains negative values.")                                         }
  OUTPUT
}
