Pkl.Hajek.U <- function(VecPk.U)
{
  if(! is.vector(VecPk.U)     ){stop("VecPk.U must be a vector.")                                                                               }
  if(any(is.na(VecPk.U))      ){stop("There are missing values in VecPk.U.")                                                                    }
  if(any(VecPk.U<=0|VecPk.U>1)){stop("There are invalid values in VecPk.U.")                                                                    }
  N                            <- length(VecPk.U)
  if(all(VecPk.U==1)          ){OUTPUT <- rep(1, times= N)                                                                                      }
  else
  {
    OUTPUT                     <- .C("Pkl_Hajek_U",
                                      as.double(VecPk.U),
                                      as.integer(N),
                                      VectMatProb = double(N*N),
                                      PACKAGE = "samplingVarEst")$VectMatProb
  }
  if(any(OUTPUT <= 0)         ){warning("Some values of 2nd order (joint) inclusion probabilities are so tiny that they are assumed zero by R.")}
  matrix(OUTPUT, nrow=N, ncol=N)
}
