Pkl.Hajek.s <- function(VecPk.s)
{
  if(! is.vector(VecPk.s)     ){stop("VecPk.s must be a vector.")                                                                               }
  if(any(is.na(VecPk.s))      ){stop("There are missing values in VecPk.s.")                                                                    }
  if(any(VecPk.s<=0|VecPk.s>1)){stop("There are invalid values in VecPk.s.")                                                                    }
  n                            <- length(VecPk.s)
  if(all(VecPk.s==1)          ){OUTPUT <- rep(1, times= n)                                                                                      }
  else
  {
    OUTPUT                     <- .C("Pkl_Hajek_s",
                                      as.double(VecPk.s),
                                      as.integer(n),
                                      VectMatProb = double(n*n),
                                      PACKAGE = "samplingVarEst")$VectMatProb
  }
  if(any(OUTPUT <= 0)         ){warning("Some values of 2nd order (joint) inclusion probabilities are so tiny that they are assumed zero by R.")}
  matrix(OUTPUT, nrow=n, ncol=n)
}
