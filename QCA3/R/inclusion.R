## This file (R/truthTable.R) is part of QCA3 package
## copyright: HUANG Ronggui 2008-2012

## implementation of inclusion algorithm for fuzzy set QCA.

suffnec.test <- function(data,outcome,conditions,type=c("suff","nec"),
                         benchmark=0.65,conf.level=0.95,...)
{
  ## outcome can be more than one
  ## only test the individual condition
  ## a separate function for the combination of conditions
  type <- match.arg(type)
  outData <- data[outcome]
  conditionsData <- data[conditions]
  ##augmentation of conditions Data by negation
  N1 <- ncol(outData)
  N2 <- ncol(conditionsData)
  ans <- matrix(nrow=N2,ncol=N1)
  rownames(ans) <- conditions
  colnames(ans) <- outcome
  if (type=="suff"){
    for (i in seq_len(N2)){
      for (j in seq_len(N1)){
        ans[i,j] <- mean(outData[j]>= conditionsData[i])
        ## add statistical test
      }
    }
  } else if (type=="nec"){
    for (i in seq_len(N2)){
      for (j in seq_len(N1)){
        ans[i,j] <- mean(outData[j] <= conditionsData[i])
        ## add statistical test
      }
    }
  }
  ans
}

##suffnec.test(KatzHauMahoney,c("EDF","EUDF","SDF","SUDF"),c("DIPF","LIEF","MTEF","SLF","SCF"),"nec")
