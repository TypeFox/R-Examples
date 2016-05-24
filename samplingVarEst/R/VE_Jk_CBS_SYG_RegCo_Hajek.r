VE.Jk.CBS.SYG.RegCo.Hajek <- function(VecY.s, VecX.s, VecPk.s, MatPkl.s)
{
  if(! is.vector(VecY.s)             ){stop("VecY.s must be a vector.")                                                     }
  if(! is.vector(VecX.s)             ){stop("VecX.s must be a vector.")                                                     }
  if(! is.vector(VecPk.s)            ){stop("VecPk.s must be a vector.")                                                    }
  if(! is.matrix(MatPkl.s)           ){stop("MatPkl.s must be a matrix.")                                                   }
  DimMat                              <- dim(MatPkl.s)
  DimMatR                             <- as.integer(DimMat[1])
  DimMatC                             <- as.integer(DimMat[2])
  if(DimMatR != DimMatC              ){stop("MatPkl.s must be a square matrix. Number of rows and columns has to be equal.")}
  n                                   <- length(VecY.s)
  if(n != length(VecPk.s)            ){stop("The lengths of VecY.s and VecPk.s are different.")                             }
  if(n != length(VecX.s)             ){stop("The lengths of VecY.s and VecX.s are different.")                              }
  if(n != DimMatR                    ){stop("The lengths of VecY.s, VecPk.s and dimensions of MatPkl.s are different.")     }
  if(any(is.na(VecPk.s))             ){stop("There are missing values in VecPk.s.")                                         }
  if(any(VecPk.s<=0|VecPk.s>1)       ){stop("There are invalid values in VecPk.s.")                                         }
  if(any(is.na(MatPkl.s))            ){stop("There are missing values in MatPkl.s.")                                        }
  if(any(MatPkl.s<=0|MatPkl.s>1)     ){stop("There are invalid values in MatPkl.s.")                                        }
  if(any(is.na(VecY.s))              ){stop("There are missing values in VecY.s.")                                          }
  if(any(is.na(VecX.s))              ){stop("There are missing values in VecX.s.")                                          }
  VecEstTheta_k                       <- .C("Est_RegCo_Hajek_Excluding_All_Elements",
                                             as.double(VecY.s),
                                             as.double(VecX.s),
                                             as.double(VecPk.s),
                                             as.integer(n),
                                             VectVarEst = double(n),
                                             PACKAGE = "samplingVarEst")$VectVarEst
  EstTheta                            <- Est.RegCo.Hajek(VecY.s, VecX.s, VecPk.s)
  Nhat                                <- .C("Est_Total_NHT",
                                             as.double(rep(1.0, times=n)),
                                             as.double(VecPk.s),
                                             as.integer(n),
                                             PointEst = double(1),
                                             PACKAGE = "samplingVarEst")$PointEst
  VecPseudo.s                         <- (1 - {1/Nhat/VecPk.s}) * (EstTheta - VecEstTheta_k)
  OUTPUT                              <- .C("VE_SYG_form",
                                             as.double(VecPseudo.s),
                                             as.double(VecPk.s),
                                             as.double(c(MatPkl.s)),
                                             as.integer(n),
                                             VarEst = double(1),
                                             PACKAGE = "samplingVarEst")$VarEst
  if(OUTPUT<0                        ){warning("The variance estimate contains negative values.")                           }
  OUTPUT
}
