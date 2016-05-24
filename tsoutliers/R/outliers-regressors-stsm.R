
##FIXME TODO SLS in outliers.regressors.stsmSS()

outliers.regressors.stsmSS <- function(pars, mo, n, weights = TRUE,
  delta = 0.7, freq = 12, n.start = 50)
{
  ioxreg <- aoxreg <- lsxreg <- tcxreg <- NULL
  
  # required to use "dind <- diff(ind)" and "seq_len(dind[i])" below
  mo <- mo[order(mo[,"ind"]),]

  # IO

##NOTE
#IO is computed as in the ARIMA case
  ind <- which(mo[,"type"] == "IO")
  lind <- length(ind)
  if (lind > 0)
  {
##NOTE
#keep this order ("w" then "ind") since "ind" is overwritten
    w <- mo[ind,"coefhat"]
    ind <- mo[ind,"ind"]
    ioxreg <- matrix(0, nrow = n, ncol = lind)
    ioxreg[n * seq.int(0, lind - 1) + ind] <- w
  }

  # common variable pi(L) I(t=1) for AO, LS and TC, named "aoxregt1"
  
  #I <- ts(rep(0, n + n.start), start = start(resid), frequency = frequency(resid))  
  #tsp(I) <- tsp(resid)
##FIXME see if KF in BSM or llmseas requires "I" to be defined as a "ts" object
  I <- rep(0, n + n.start)
  I[n.start] <- 1
  tmp <- KFKSDS::KF(I, pars)
##NOTE
#aoxregt1 is "f" in tstatsSTSM()
#pi(L) I(t=1)
  aoxregt1 <- tmp$v[-seq.int(n.start-1)]
##FIXME see how to avoid remove last element this way
  aoxregt1 <- aoxregt1[-length(aoxregt1)]

  # AO

  ind <- which(mo[,"type"] == "AO")
  lind <- length(ind)
  if (lind > 0)
  {
##NOTE
#keep this order ("w" then "ind") since "ind" is overwritten
    w <- mo[ind,"coefhat"]
    ind <- mo[ind,"ind"]
    dind <- diff(ind)
    
    aoxreg <- matrix(0, nrow = n, ncol = lind)
    id <- seq.int(ind[1], n)    
    
    for (i in seq_along(ind))
    {
      #id <- seq.int(ind[i], n)
      #d[id,i] <- b[seq_along(id)]
      #id <- id[-seq_len(ind[i])]
      aoxreg[id,i] <- w[i] * aoxregt1[seq_along(id)]
      if (i != lind)
       id <- id[-seq_len(dind[i])]
    }
  }
  
  # LS
  
  ind <- which(mo[,"type"] == "LS")
  lind <- length(ind)
  if (lind > 0)
  {
##NOTE
#keep this order ("w" then "ind") since "ind" is overwritten
    w <- mo[ind,"coefhat"]
    ind <- mo[ind,"ind"]
    dind <- diff(ind)
    
    lsxreg <- matrix(0, nrow = n, ncol = lind)
    id <- seq.int(ind[1], n)    
    
    for (i in seq_along(ind))
    {
      #id <- seq.int(ind[i], n)
      #d[id,i] <- b[seq_along(id)]
      #id <- id[-seq_len(ind[i])]
      lsxreg[id,i] <- w[i] * diffinv(aoxregt1[seq_along(id)])[-1]
      if (i != lind)
       id <- id[-seq_len(dind[i])]
    }
  }

  # TC
  
  ind <- which(mo[,"type"] == "TC")
  lind <- length(ind)
  if (lind > 0)
  {
##NOTE
#keep this order ("w" then "ind") since "ind" is overwritten
    w <- mo[ind,"coefhat"]
    ind <- mo[ind,"ind"]
    dind <- diff(ind)
    
    tcxreg <- matrix(0, nrow = n, ncol = lind)
    id <- seq.int(ind[1], n)    
    
    for (i in seq_along(ind))
    {
      #id <- seq.int(ind[i], n)
      #d[id,i] <- b[seq_along(id)]
      #id <- id[-seq_len(ind[i])]
      tcxreg[id,i] <- w[i] * diffinv(aoxregt1[seq_along(id)])[-1]
      if (i != lind)
       id <- id[-seq_len(dind[i])]
    }
  }

##FIXME see give colnames for clarity but if used just to sum rows it may not worth doing it
##NOTE the columns are not in the same order as the rows in "mo" (be careful if weights=FALSE
#is used and then the matrix is used somehow, do not assume same order as "mo")
  cbind(ioxreg, aoxreg, lsxreg, tcxreg)
}
