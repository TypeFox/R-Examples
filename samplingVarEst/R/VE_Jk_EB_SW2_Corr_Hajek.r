VE.Jk.EB.SW2.Corr.Hajek <- function(VecY.s, VecX.s, VecPk.s, nII, VecPi.s, VecCluLab.s, VecCluSize.s)
{
  if(! is.vector(VecY.s)                    ){stop("VecY.s must be a vector.")                                                                                               }
  if(! is.vector(VecX.s)                    ){stop("VecX.s must be a vector.")                                                                                               }
  if(! is.vector(VecPk.s)                   ){stop("VecPk.s must be a vector.")                                                                                              }
  if(! is.vector(VecPi.s)                   ){stop("VecPi.s must be a vector.")                                                                                              }
  if(! is.vector(VecCluLab.s)               ){stop("VecCluLab.s must be a vector.")                                                                                          }
  if(! is.vector(VecCluSize.s)              ){stop("VecCluSize.s must be a vector.")                                                                                         }
  if(! is.integer(VecCluSize.s)             ){stop("VecCluSize.s must be a vector of integers.")                                                                             }
  if(! is.vector(VecCluLab.s)               ){stop("VecCluLab.s must be a vector.")                                                                                          }
  if(! is.integer(VecCluLab.s)              ){stop("VecCluLab.s must be a vector of integers.")                                                                              }
  if(any(is.na(VecPk.s))                    ){stop("There are missing values in VecPk.s.")                                                                                   }
  if(any(VecPk.s<=0|VecPk.s>1)              ){stop("There are invalid values in VecPk.s.")                                                                                   }
  if(any(is.na(VecPi.s))                    ){stop("There are missing values in VecPi.s.")                                                                                   }
  if(any(VecPi.s<=0|VecPi.s>1)              ){stop("There are invalid values in VecPi.s.")                                                                                   }
  if(any(is.na(VecY.s))                     ){stop("There are missing values in VecY.s.")                                                                                    }
  if(any(is.na(VecX.s))                     ){stop("There are missing values in VecX.s.")                                                                                    }
  if(length(nII) != 1                       ){stop(paste("nII must be a integer scalar, i.e. a vector of length 1 valued from 1 to", as.character(.Machine$integer.max),"."))}
  n                                          <- length(VecY.s)
  VecPi.sI                                   <- VecPi.s[!duplicated(VecCluLab.s)]
  VecCluSize.sI                              <- VecCluSize.s[!duplicated(VecCluLab.s)]
  nI                                         <- length(VecPi.sI)
  DoublenI                                   <- as.double(nI)
  if(n != length(VecPk.s)                   ){stop("The lengths of VecY.s and VecPk.s are different.")                                                                       }
  if(n != length(VecX.s)                    ){stop("The lengths of VecY.s and VecX.s are different.")                                                                        }
  if(n != length(VecCluLab.s)               ){stop("The lengths of VecY.s and VecCluLab.s are different.")                                                                   }
  if(n != as.integer(nI*nII)                ){stop("The lengths of VecY.s and (nI * nII) are different.")                                                                    }
  if(any(VecCluSize.s < nII)                ){stop("The value of nII is larger than the minimum cluster size in VecCluSize.s.")                                              }
  EstTheta                                   <- Est.Corr.Hajek(VecY.s, VecX.s, VecPk.s)
  VarianceClusters                           <- .C("VE_Jk_EB_SW2_Corr_Hajek_Clu",
                                                    as.double(EstTheta),
                                                    as.double(VecY.s),
                                                    as.double(VecX.s),
                                                    as.double(VecPk.s),
                                                    as.integer(n),
                                                    as.integer(nI),
                                                    as.double(DoublenI),
                                                    as.integer(VecCluLab.s),
                                                    as.integer(VecCluLab.s[!duplicated(VecCluLab.s)]),
                                                    as.double(1 - VecPi.sI),
                                                    as.double(1 - VecPi.sI * nII/(nII-1)*(VecCluSize.sI-1)/VecCluSize.sI ),
                                                    Est_d = double(1),
                                                    VarianceClusters = double(1),
                                                    PACKAGE = "samplingVarEst")$VarianceClusters
  VecEstTheta_k                              <- .C("Est_Corr_Hajek_Excluding_All_Elements",
                                                    as.double(VecY.s),
                                                    as.double(VecX.s),
                                                    as.double(VecPk.s),
                                                    as.integer(n),
                                                    VectVarEst = double(n),
                                                    PACKAGE = "samplingVarEst")$VectVarEst
  SQVecPseudo.s                              <- ( (n-1)/n * (EstTheta - VecEstTheta_k) )^2
  VectClusterSumSqPseudoValues               <- rowsum(SQVecPseudo.s, group = VecCluLab.s)
  VectPhiClu.sI                              <- VecPi.sI * nII/(nII -1) * (1 - nII/VecCluSize.sI)
  VarianceObservations                       <- sum(VectPhiClu.sI * VectClusterSumSqPseudoValues)
  OUTPUT                                     <- VarianceClusters + VarianceObservations
  if(OUTPUT<0                               ){warning("The variance estimate contains negative values.")                                                                     }
  OUTPUT
}
