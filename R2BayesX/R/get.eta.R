get.eta <-
function(data)
{
  etaspec <- c("eta","linpred","pmean_pred","pqu2p5_pred","pqu10_pred",
    "pmed_pred","pqu90_pred","pqu97p5_pred","pmean_mu","pqu2p5_mu",
    "pqu10_mu","pmed_mu","pqu90_mu","pqu97p5_mu","mu")
  nd <- names(data)
  eta <- enam <- NULL
  for(e in etaspec) 
    for(n in nd)
      if(!is.na(which <- pmatch(e, n))) {
        eta <- cbind(eta, as.matrix(data[n], mode = "numeric")) 
        enam <- c(enam, n)
      }
  if(!is.null(eta)) {
    colnames(eta) <- enam
    rownames(eta) <- 1L:nrow(eta)
    storage.mode(eta) <- "numeric"
  }
  if(!is.null(eta))
    eta <- na.omit(eta)

  return(eta)
}

