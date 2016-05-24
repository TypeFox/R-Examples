Pk.PropNorm.U <- function(n, VecMOS.U)
{
  if(! is.vector(VecMOS.U)){stop("VecMOS.U must be a vector.")                                                                                     }
  if(any(is.na(VecMOS.U)) ){stop("There are missing values in VecMOS.U.")                                                                          }
  if(any(VecMOS.U<=0)     ){warning("There are zero or negative values in VecMOS.U.")                                                              }
  if(n%%1 != 0            ){stop("n must be an integer or a double-precision scalar with zero-valued fractional part.")                            }
  N                        <- length(VecMOS.U)
  OUTPUT                   <- .C("Pk_PropNorm_U",
                                      as.integer(n),
                                      as.integer(N),
                                      as.double(VecMOS.U),
                                      VecPk = double(N),
                                      PACKAGE = "samplingVarEst")$VecPk
  if(any(OUTPUT <= 0)     ){warning("Some normalised values of the 1st order inclusion probabilities are so tiny that they are assumed zero by R.")}
  as.vector(OUTPUT)
}
