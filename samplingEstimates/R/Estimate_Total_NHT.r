Estimate.Total.NHT <- function(MatY.s, VecWk.s, VarEst= "SYG", MatPkl.s= NULL, PopSize= NULL, VecStratLb.s= NULL, VecStratSize.s= NULL, ShowStrata= FALSE, VecDomainLb.s= NULL)
{
  ###Checking the input and displaying messages

  #Checking MatY.s
  if(  is.vector(MatY.s)             ){tempname <- as.character(deparse(substitute(MatY.s))); MatY.s <- as.matrix(MatY.s, ncol =1); colnames(MatY.s) <- tempname                 }
  if(any(is.na(MatY.s))              ){stop("There are missing values in MatY.s.")                                                                                               }
  n                                   <- dim(MatY.s)[1]
  Q                                   <- dim(MatY.s)[2]

  #Checking VecWk.s
  if(! is.vector(VecWk.s)            ){stop("VecWk.s must be a vector.")                                                                                                         }
  if(any(is.na(VecWk.s))             ){stop("There are missing values in VecWk.s.")                                                                                              }
  if(any(VecWk.s<1)                  ){stop("There are invalid values in VecWk.s, some entries are < 1.")                                                                        }
  if(n != length(VecWk.s)            ){stop("Lengths of columns of MatY.s and the length of VecWk.s are different.")                                                             }

  #Computing some quantities
  Nhat                                <- sum(VecWk.s)
  fhat                                <- n/Nhat
  VecPk.s                             <- 1.0/VecWk.s

  #Checking MatPkl.s
  if(!is.null(MatPkl.s))if(any(is.na(MatPkl.s))){stop("There are missing values in MatPkl.s.")                                                                                   }
  if(!is.null(MatPkl.s))if(!is.matrix(MatPkl.s)){stop("MatPkl.s must be a matrix.")                                                                                              }
  if(is.null(MatPkl.s)               ){VarEst <- "Hajek"; cat("Note: MatPkl.s, the joint inclusion probabilities are not provided, VarEst is set to \"Hajek\".\n")               }
  if(     VarEst == "Hajek"          ){cat("Note: when using Hajek's approximations, a high-entropy sampling design is assumed.\n")                                              }

  #Checking PopSize
  if(!is.null(PopSize)               ){f <- n/PopSize                                                                                                                            }

  #First steps of preparing the output
  VecVarName                          <- colnames(MatY.s)
  OUTPUT                              <- data.frame(cbind(Statistic = rep("Total.NHT", times= Q), VariableName = VecVarName))

  #Checking if stratification is present
  if(is.null(VecStratLb.s) & !is.null(VecStratSize.s)){stop("The argument VecStratSize.s is provided but the argument VecStratLb.s is missing.")                                 }
  if(!is.null(VecStratLb.s) & is.null(VecStratSize.s)){stop("The argument VecStratLb.s is provided but the argument VecStratSize.s is missing.")                                 }
  FlagNoStratification                <- is.null(VecStratLb.s) & is.null(VecStratSize.s)

  if(FlagNoStratification) #Sub-routine if there is no stratification
  {
    for(q in (1:Q))
    {
      OUTPUT$Estimate[q]              <- Est.Total.NHT(MatY.s[,q], VecPk.s)
      if(     VarEst == "HT"         ){OUTPUT$Variance[q] <- VE.HT.Total.NHT(   MatY.s[,q], VecPk.s, MatPkl.s)                                                                   }
      else if(VarEst == "SYG"        ){OUTPUT$Variance[q] <- VE.SYG.Total.NHT(  MatY.s[,q], VecPk.s, MatPkl.s)                                                                   }
      else if(VarEst == "Hajek"      ){OUTPUT$Variance[q] <- VE.Hajek.Total.NHT(MatY.s[,q], VecPk.s          )                                                                   }
      else                            {stop("The argument VarEst must be: \"HT\", \"SYG\" or \"Hajek\". If omitted, default is \"SYG\".")                                        }
      OUTPUT$StdErr[q]                <- sqrt(OUTPUT$Variance[q])
      OUTPUT$AbsErr[q]                <- OUTPUT$StdErr[q] * 1.959
      OUTPUT$LInfCI95[q]              <- OUTPUT$Estimate[q] - OUTPUT$AbsErr[q]; if(OUTPUT$LInfCI95[q]<0){OUTPUT$LInfCI95[q] <- 0; cat("Note: LInfCI95 < 0, it is set to zero.\n")}
      OUTPUT$LSupCI95[q]              <- OUTPUT$Estimate[q] + OUTPUT$AbsErr[q]
      OUTPUT$Range95[q]               <- OUTPUT$LSupCI95[q] - OUTPUT$LInfCI95[q]
      OUTPUT$PctCVE[q]                <- round(100.0*OUTPUT$StdErr[q]/OUTPUT$Estimate[q], digits= 3)
      if(!is.null(PopSize)           ){OUTPUT$DEff[q] <- round(OUTPUT$Variance[q]/VE.SYG.Total.NHT(MatY.s[,q],rep(f,   times=n),Pkl.Hajek.s(rep(f,   times=n))), digits = 5)     }
      else                            {OUTPUT$DEff[q] <- round(OUTPUT$Variance[q]/VE.SYG.Total.NHT(MatY.s[,q],rep(fhat,times=n),Pkl.Hajek.s(rep(fhat,times=n))), digits = 5)     }
    }
    OUTPUT$n                          <- n
    OUTPUT$Nhat                       <- round(Nhat, digits= 2)
    OUTPUT$fhat                       <- round(fhat, digits= 4)
    if(!is.null(PopSize)             ){OUTPUT$N <- PopSize; OUTPUT$f <- round(f, digits= 4)                                                                                      }
  }
  else if(!FlagNoStratification) #Sub-routine if there is stratification
  {
    if(any(is.na(VecStratLb.s))      ){stop("There are missing values in VecStratLb.s.")                                                                                         }
    if(any(is.na(VecStratSize.s))    ){stop("There are missing values in VecStratSize.s.")                                                                                       }
    if(n != length(VecStratLb.s)     ){stop("Lengths of columns of MatY.s and the length of VecStratLb.s are different.")                                                        }
    if(n != length(VecStratSize.s)   ){stop("Lengths of columns of MatY.s and the length of VecStratSize.s are different.")                                                      }
    if(is.factor(VecStratLb.s)       ){VecStratLb.s <- as.character(VecStratLb.s)                                                                                                }
    VecStrataLbls.H                   <- unique(VecStratLb.s)
    H                                 <- length(VecStrataLbls.H)
    for(q in (1:Q))
    {
      OUTPUTSTRATA                    <- data.frame(cbind(Statistic = rep("Total.NHT", times= H), VariableName = VecVarName[q]))
      for(h in (1:H))
      {
        OUTPUTSTRATA$h[h]             <- h
        OUTPUTSTRATA$Stratum[h]       <- VecStrataLbls.H[h]
        OUTPUTSTRATA$Estimate[h]      <- Est.Total.NHT(MatY.s[VecStratLb.s == VecStrataLbls.H[h], q], VecPk.s[VecStratLb.s == VecStrataLbls.H[h]])
        if(     VarEst== "HT"        ){OUTPUTSTRATA$Variance[h] <- VE.HT.Total.NHT(   MatY.s[VecStratLb.s== VecStrataLbls.H[h],q], VecPk.s[VecStratLb.s== VecStrataLbls.H[h]],MatPkl.s[VecStratLb.s== VecStrataLbls.H[h],VecStratLb.s== VecStrataLbls.H[h]])}
        else if(VarEst== "SYG"       ){OUTPUTSTRATA$Variance[h] <- VE.SYG.Total.NHT(  MatY.s[VecStratLb.s== VecStrataLbls.H[h],q], VecPk.s[VecStratLb.s== VecStrataLbls.H[h]],MatPkl.s[VecStratLb.s== VecStrataLbls.H[h],VecStratLb.s== VecStrataLbls.H[h]])}
        else if(VarEst== "Hajek"     ){OUTPUTSTRATA$Variance[h] <- VE.Hajek.Total.NHT(MatY.s[VecStratLb.s== VecStrataLbls.H[h],q], VecPk.s[VecStratLb.s== VecStrataLbls.H[h]])   }
        else                          {stop("The argument VarEst must be: \"HT\", \"SYG\" or \"Hajek\". If omitted, default is \"SYG\".")                                        }
        OUTPUTSTRATA$StdErr[h]        <- sqrt(OUTPUTSTRATA$Variance[h])
        OUTPUTSTRATA$AbsErr[h]        <- OUTPUTSTRATA$StdErr[h] * 1.959
        OUTPUTSTRATA$LInfCI95[h]      <- OUTPUTSTRATA$Estimate[h] - OUTPUTSTRATA$AbsErr[h]; if(OUTPUTSTRATA$LInfCI95[h]<0){OUTPUTSTRATA$LInfCI95[h] <- 0                         }
        OUTPUTSTRATA$LSupCI95[h]      <- OUTPUTSTRATA$Estimate[h] + OUTPUTSTRATA$AbsErr[h]
        OUTPUTSTRATA$Range95[h]       <- OUTPUTSTRATA$LSupCI95[h] - OUTPUTSTRATA$LInfCI95[h]
        OUTPUTSTRATA$PctCVE[h]        <- round(100.0*OUTPUTSTRATA$StdErr[h]/OUTPUTSTRATA$Estimate[h], digits= 3)
        MaxStratSizeWithinh           <- max(VecStratSize.s[VecStratLb.s== VecStrataLbls.H[h]])
        MinStratSizeWithinh           <- min(VecStratSize.s[VecStratLb.s== VecStrataLbls.H[h]])
        if( MaxStratSizeWithinh != MinStratSizeWithinh){stop("Something is wrong with VecStratSize.s. Strata sizes for observations within the same stratum are different.")     }
        Nh                            <- MaxStratSizeWithinh
        nh                            <- sum(VecStratLb.s== VecStrataLbls.H[h])
        fh                            <- nh/Nh
        OUTPUTSTRATA$DEff[h]          <- round(OUTPUTSTRATA$Variance[h]/VE.SYG.Total.NHT(MatY.s[VecStratLb.s== VecStrataLbls.H[h],q],rep(fh,times=nh),Pkl.Hajek.s(rep(fh,times=nh))), digits= 5)
        OUTPUTSTRATA$nh[h]            <- nh
        OUTPUTSTRATA$Nh[h]            <- Nh
        OUTPUTSTRATA$fh[h]            <- round(fh, digits= 5)
      }
      if(!is.null(PopSize))if(PopSize != sum(OUTPUTSTRATA$Nh)){stop("Something is wrong with VecStratSize.s. The sum of the identified strata sizes does not equal PopSize.")    }
      OUTPUTSTRATA$Wh                 <- OUTPUTSTRATA$Nh/sum(OUTPUTSTRATA$Nh)
      if(ShowStrata                  ){cat("\n"); print(OUTPUTSTRATA, row.names = FALSE)                                                                                         }
      OUTPUT$Estimate[q]              <- sum(OUTPUTSTRATA$Estimate)
      OUTPUT$Variance[q]              <- sum(OUTPUTSTRATA$Variance)
      OUTPUT$StdErr[q]                <- sqrt(OUTPUT$Variance[q])
      OUTPUT$AbsErr[q]                <- OUTPUT$StdErr[q] * 1.959
      OUTPUT$LInfCI95[q]              <- OUTPUT$Estimate[q] - OUTPUT$AbsErr[q]; if(OUTPUT$LInfCI95[q]<0){OUTPUT$LInfCI95[q] <- 0; cat("Note: LInfCI95 < 0, it is set to zero.\n")}
      OUTPUT$LSupCI95[q]              <- OUTPUT$Estimate[q] + OUTPUT$AbsErr[q]
      OUTPUT$Range95[q]               <- OUTPUT$LSupCI95[q] - OUTPUT$LInfCI95[q]
      OUTPUT$PctCVE[q]                <- round(100.0*OUTPUT$StdErr[q]/OUTPUT$Estimate[q], digits = 3)
      if(!is.null(PopSize)           ){OUTPUT$DEff[q] <- round(OUTPUT$Variance[q]/VE.SYG.Total.NHT(MatY.s[,q],rep(f,   times=n),Pkl.Hajek.s(rep(f,   times=n))), digits = 5)     }
      else                            {OUTPUT$DEff[q] <- round(OUTPUT$Variance[q]/VE.SYG.Total.NHT(MatY.s[,q],rep(fhat,times=n),Pkl.Hajek.s(rep(fhat,times=n))), digits = 5)     }
    }
    OUTPUT$n                          <- n
    OUTPUT$Nhat                       <- round(Nhat, digits= 2)
    OUTPUT$fhat                       <- round(fhat, digits= 5)
    if(!is.null(PopSize)             ){OUTPUT$N <- PopSize; OUTPUT$f <- round(f, digits= 5)                                                                                      }
  }

  #Checking if domains' results are requested for on-screen display
  if(! is.null(VecDomainLb.s)        )
  {
    if(n != length(VecDomainLb.s)    ){stop("Lengths of columns of MatY.s and the length of VecDomainLb.s are different.")                                                       }
    if(any(is.na(VecDomainLb.s))     ){stop("There are missing values in VecDomainLb.s.")                                                                                        }
    if(is.factor(VecDomainLb.s)      ){VecDomainLb.s <- as.character(VecDomainLb.s)                                                                                              }
    VecDomainLbls.D                   <- unique(VecDomainLb.s)
    D                                 <- length(VecDomainLbls.D)
    for(q in (1:Q))
    {
      OUTPUTDOMAIN                    <- data.frame(cbind(Statistic = rep("Total.NHT", times= D), VariableName = VecVarName[q]))
      for(d in (1:D))
      {
        OUTPUTDOMAIN$d[d]             <- d
        OUTPUTDOMAIN$Domain[d]        <- VecDomainLbls.D[d]
        OUTPUTDOMAIN$Estimate[d]      <- Est.Total.NHT(MatY.s[VecDomainLb.s == VecDomainLbls.D[d], q], VecPk.s[VecDomainLb.s == VecDomainLbls.D[d]])
        if(     VarEst== "HT"        ){OUTPUTDOMAIN$Variance[d] <- VE.HT.Total.NHT(   MatY.s[VecDomainLb.s== VecDomainLbls.D[d],q], VecPk.s[VecDomainLb.s== VecDomainLbls.D[d]],MatPkl.s[VecDomainLb.s== VecDomainLbls.D[d],VecDomainLb.s== VecDomainLbls.D[d]])}
        else if(VarEst== "SYG"       ){OUTPUTDOMAIN$Variance[d] <- VE.SYG.Total.NHT(  MatY.s[VecDomainLb.s== VecDomainLbls.D[d],q], VecPk.s[VecDomainLb.s== VecDomainLbls.D[d]],MatPkl.s[VecDomainLb.s== VecDomainLbls.D[d],VecDomainLb.s== VecDomainLbls.D[d]])}
        else if(VarEst== "Hajek"     ){OUTPUTDOMAIN$Variance[d] <- VE.Hajek.Total.NHT(MatY.s[VecDomainLb.s== VecDomainLbls.D[d],q], VecPk.s[VecDomainLb.s== VecDomainLbls.D[d]]) }
        else                          {stop("The argument VarEst must be: \"HT\", \"SYG\" or \"Hajek\". If omitted, default is \"SYG\".")                                        }
        OUTPUTDOMAIN$StdErr[d]        <- sqrt(OUTPUTDOMAIN$Variance[d])
        OUTPUTDOMAIN$AbsErr[d]        <- OUTPUTDOMAIN$StdErr[d] * 1.959
        OUTPUTDOMAIN$LInfCI95[d]      <- OUTPUTDOMAIN$Estimate[d] - OUTPUTDOMAIN$AbsErr[d]; if(OUTPUTDOMAIN$LInfCI95[d]<0){OUTPUTDOMAIN$LInfCI95[d] <- 0                         }
        OUTPUTDOMAIN$LSupCI95[d]      <- OUTPUTDOMAIN$Estimate[d] + OUTPUTDOMAIN$AbsErr[d]
        OUTPUTDOMAIN$Range95[d]       <- OUTPUTDOMAIN$LSupCI95[d] - OUTPUTDOMAIN$LInfCI95[d]
        OUTPUTDOMAIN$PctCVE[d]        <- round(100.0*OUTPUTDOMAIN$StdErr[d]/OUTPUTDOMAIN$Estimate[d], digits= 3)
        nd                            <- sum(VecDomainLb.s== VecDomainLbls.D[d])
        Ndhat                         <- sum(VecWk.s[VecDomainLb.s== VecDomainLbls.D[d]])
        fdhat                         <- nd/Ndhat
        OUTPUTDOMAIN$DEff[d]          <- round(OUTPUTDOMAIN$Variance[d]/VE.SYG.Total.NHT(MatY.s[VecDomainLb.s== VecDomainLbls.D[d],q],rep(fdhat,times=nd),Pkl.Hajek.s(rep(fdhat,times=nd))), digits= 5)
        OUTPUTDOMAIN$nd[d]            <- nd
        OUTPUTDOMAIN$Ndhat[d]         <- Ndhat
        OUTPUTDOMAIN$fdhat[d]         <- round(fdhat, digits= 5)
        if(!is.null(PopSize)         ){OUTPUTDOMAIN$Wdhat[d] <- Ndhat/PopSize                                                                                                    }
        else                          {OUTPUTDOMAIN$Wdhat[d] <- Ndhat/Nhat                                                                                                       }
      }
      cat("\n"); print(OUTPUTDOMAIN, row.names = FALSE)
    }
  }

  #Returning the main result
  cat("\n")
  return(OUTPUT)
}
