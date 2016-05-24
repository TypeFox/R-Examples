Projector <- function() {
  Rcmdr <- options()$Rcmdr
  Projector.options <-
    list(log.font.size = 15,
         log.width = 54,
         log.height = 6, 
         output.height = 18,
         scale.factor = 1.4)
  Rcmdr[names(Projector.options)] <- Projector.options
  options(Rcmdr=Rcmdr)
  trellis.par.set(list(superpose.symbol=list(
                         pch=rep(16,
                           length(trellis.par.get("superpose.symbol")$pch))),
                       plot.symbol=list(pch=16)))
  par(pch=16)
  
  putRcmdr("autoRestart", TRUE)
  closeCommander(ask=FALSE)
  Commander()
}

## source("~/HH-R.package/RcmdrPlugin.HH/R/Projector.R")
