VE.EB.HT.Ratio <- function(VecY.s, VecX.s, VecPk.s, MatPkl.s, VecAlpha.s = rep(1, times= length(VecPk.s)))
{
  if(! is.vector(VecY.s)                ){stop("VecY.s must be a vector.")                                                        }
  if(! is.vector(VecX.s)                ){stop("VecX.s must be a vector.")                                                        }
  if(! is.vector(VecPk.s)               ){stop("VecPk.s must be a vector.")                                                       }
  if(! is.vector(VecAlpha.s)            ){stop("VecAlpha.s must be a vector.")                                                    }
  if(! is.matrix(MatPkl.s)              ){stop("MatPkl.s must be a matrix.")                                                      }
  DimMat                                 <- dim(MatPkl.s)
  DimMatR                                <- as.integer(DimMat[1])
  DimMatC                                <- as.integer(DimMat[2])
  if(DimMatR != DimMatC                 ){stop("MatPkl.s must be a square matrix. Number of rows and columns has to be equal.")   }
  n                                      <- length(VecY.s)
  if(n != length(VecPk.s)               ){stop("The lengths of VecY.s and VecPk.s are different.")                                }
  if(n != length(VecX.s)                ){stop("The lengths of VecY.s and VecX.s are different.")                                 }
  if(n != length(VecAlpha.s)            ){stop("The lengths of VecY.s and VecAlpha.s are different.")                             }
  if(n != DimMatR                       ){stop("The lengths of VecY.s, VecPk.s and dimensions of MatPkl.s are different.")        }
  if(any(is.na(VecPk.s))                ){stop("There are missing values in VecPk.s.")                                            }
  if(any(is.na(VecAlpha.s))             ){stop("There are missing values in VecAlpha.s.")                                         }
  if(any(VecPk.s<=0|VecPk.s>1)          ){stop("There are invalid values in VecPk.s.")                                            }
  if(any(VecAlpha.s<0)                  ){stop("There are invalid values in VecAlpha.s.")                                         }
  if(any(VecAlpha.s==0)                 ){warning("There are zero values in VecAlpha.s. It is valid, but estimates get unstable.")}
  if(any(is.na(MatPkl.s))               ){stop("There are missing values in MatPkl.s.")                                           }
  if(any(MatPkl.s<=0|MatPkl.s>1)        ){stop("There are invalid values in MatPkl.s.")                                           }
  if(any(is.na(VecY.s))                 ){stop("There are missing values in VecY.s.")                                             }
  if(any(is.na(VecX.s))                 ){stop("There are missing values in VecX.s.")                                             }
  if(any(VecX.s==0)                     ){warning("There are zero values in the denominator variable VecX.s.")                    }
  VecRhok.s                              <- (1/VecPk.s)^(1-VecAlpha.s)
  VecEstTheta_k                          <- .C("Est_Ratio_Excluding_All_ElementsByRho",
                                                as.double(VecY.s),
                                                as.double(VecX.s),
                                                as.double(VecPk.s),
                                                as.double(VecRhok.s),
                                                as.integer(n),
                                                VectVarEst = double(n),
                                                PACKAGE = "samplingVarEst")$VectVarEst
  EstTheta                               <- Est.Ratio(VecY.s, VecX.s, VecPk.s)
  VecPseudo.s                            <- (EstTheta - VecEstTheta_k)/(VecPk.s*VecRhok.s)
  OUTPUT                                 <- .C("VE_HT_form",
                                                as.double(VecPseudo.s),
                                                as.double(VecPk.s),
                                                as.double(c(MatPkl.s)),
                                                as.integer(n),
                                                VarEst = double(1),
                                                PACKAGE = "samplingVarEst")$VarEst
  if(OUTPUT<0                           ){warning("The variance estimate contains negative values.")                              }
  OUTPUT
}
