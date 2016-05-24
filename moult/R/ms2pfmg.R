ms2pfmg <- function (ms, fm) 
{ pfmg <- numeric()
  ps <- c(0, 0.125, 0.375, 0.625, 0.875, 1)
  totalmass <- sum(fm)
  j <- length(fm)

  if(is.numeric(ms))
    scores <- format(ms,scientific=FALSE,trim=TRUE)
  else
    scores <- ms
  for (n in 1:length(ms))
    { pm <- 0
      for (i in 1:j)
        { pm <- pm + ps[as.numeric(substr(scores[n], i, i)) + 1] * fm[i]/totalmass
        }
      pfmg <- c(pfmg, pm)
    }
  return(pfmg)
}
