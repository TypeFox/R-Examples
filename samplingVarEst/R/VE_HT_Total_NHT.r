VE.HT.Total.NHT <- function(VecY.s, VecPk.s, MatPkl.s)
{
  if(! is.vector(VecY.s)             ){stop("VecY.s must be a vector.")                                                     }
  if(! is.vector(VecPk.s)            ){stop("VecPk.s must be a vector.")                                                    }
  if(! is.matrix(MatPkl.s)           ){stop("MatPkl.s must be a matrix.")                                                   }
  DimMat                              <- dim(MatPkl.s)
  DimMatR                             <- as.integer(DimMat[1])
  DimMatC                             <- as.integer(DimMat[2])
  if(DimMatR != DimMatC              ){stop("MatPkl.s must be a square matrix. Number of rows and columns has to be equal.")}
  n                                   <- length(VecY.s)
  if(n != length(VecPk.s)            ){stop("The lengths of VecY.s and VecPk.s are different.")                             }
  if(n != DimMatR                    ){stop("The lengths of VecY.s, VecPk.s and dimensions of MatPkl.s are different.")     }
  if(any(is.na(VecPk.s))             ){stop("There are missing values in VecPk.s.")                                         }
  if(any(VecPk.s<=0|VecPk.s>1)       ){stop("There are invalid values in VecPk.s.")                                         }
  if(any(is.na(MatPkl.s))            ){stop("There are missing values in MatPkl.s.")                                        }
  if(any(MatPkl.s<=0|MatPkl.s>1)     ){stop("There are invalid values in MatPkl.s.")                                        }
  if(any(is.na(VecY.s))              ){stop("There are missing values in VecY.s.")                                          }
  OUTPUT                              <- .C("VE_HT_Total_NHT",
                                             as.double(VecY.s),
                                             as.double(VecPk.s),
                                             as.double(c(MatPkl.s)),
                                             as.integer(n),
                                             VarEst = double(1),
                                             PACKAGE = "samplingVarEst")$VarEst
  if(OUTPUT<0                        ){warning("The variance estimate contains negative values.")                           }
  OUTPUT
}
