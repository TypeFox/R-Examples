gammaHRF <- function(TR, paras=NULL, len.seconds=32, onset.seconds=0) {
  if (is.null(paras)) paras <- c(6,16,1,1,6)
  dt <- TR/16
  u <- 0:(len.seconds/dt) - onset.seconds/dt
  hrf <- dgamma(u, paras[1]/paras[3], dt/paras[3]) - dgamma(u, paras[2]/paras[4], dt/paras[4])/paras[5]
  hrf <- hrf[(0:(len.seconds/TR))*16+1]
  hrf <- hrf/sum(hrf)
  return(hrf)
}
