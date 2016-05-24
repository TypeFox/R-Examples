cutLMms <- function(msdata, LM, cens)
{
  tmat <- attr(msdata, "trans")
  msdata$Tentry <- msdata$Tstart
  msdata <- msdata[msdata$Tstop>LM,]
  msdata$Tstart[msdata$Tstart<LM] <- LM
  if (!missing(cens)) {
    msdata$status[msdata$Tstop>cens] <- 0
    msdata$Tstop[msdata$Tstop>cens] <- cens
    msdata <- msdata[msdata$Tstop>=msdata$Tstart,]
  }
  msdata$time <- msdata$Tstop - msdata$Tstart
  attr(msdata, "trans") <- tmat
  return(msdata)
}
