dftest <-
function(y){
  out<-summary(ur.df(y, lags=0, type="drift"))
  DFStat <- c(slot(out, "teststat")[1,1])
  names(DFStat) <- "Dickey-Fuller"
  pv <- slot(out, "cval")[1,]
  list(tau =DFStat , criticalValues=pv)
}

