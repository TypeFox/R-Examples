Est.Total.NHT <- function(VecY.s, VecPk.s)
{
  if(! is.vector(VecY.s)             ){stop("VecY.s must be a vector.")                        }
  if(! is.vector(VecPk.s)            ){stop("VecPk.s must be a vector.")                       }
  if(any(is.na(VecPk.s))             ){stop("There are missing values in VecPk.s.")            }
  if(any(VecPk.s<=0|VecPk.s>1)       ){stop("There are invalid values in VecPk.s.")            }
  if(any(is.na(VecY.s))              ){stop("There are missing values in VecY.s.")             }
  n                                   <- length(VecY.s)
  if(n != as.integer(length(VecPk.s))){stop("The lengths of VecY.s and VecPk.s are different.")}
  .C("Est_Total_NHT", 
      as.double(VecY.s), 
      as.double(VecPk.s), 
      as.integer(n), 
      PointEst = double(1), 
      PACKAGE = "samplingVarEst")$PointEst
}
