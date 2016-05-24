xsect <- function(msdata, xtime=0)
{
  msd <- msdata[msdata$Tstart<=xtime & msdata$Tstop>xtime,]
  msd <- msd[order(msd$id, msd$trans, msd$Tstart),]
  msd <- msd[!duplicated(msd$id),]
  idstate <- data.frame(id=msd$id, state=msd$from)
  tmat <- attr(msdata, "trans")
  K <- nrow(tmat)
  msd$from <- factor(msd$from, levels=1:K, labels=1:K)
  tbl <- table(msd$from)
  atrisk <- sum(tbl)
  prop <- tbl/atrisk
  res <- idstate
  attr(res, "atrisk") <- tbl
  attr(res, "prop") <- prop
  return(res)
}
