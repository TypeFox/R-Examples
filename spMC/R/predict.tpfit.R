predict.tpfit <-
function(object, lags, ...) {
  # Transition probabilities prediction for spatial Markov Chain in 1D
  #
  #    lags vector of lags
  #  object list with a matrix of estimated transition rates

  storage.mode(lags) <- "double"
  mydim <- c(dim(object$coefficients), length(lags))
  mypred <- array(0, dim = mydim)
  mypred <- .C('predTPFIT', coefficients = as.double(object$coefficients), 
               prop = as.double(object$prop), lags = as.double(lags), 
               mydim = as.integer(mydim), mypred = as.double(mypred), 
               PACKAGE = "spMC")$mypred

  res <- list()
  res$Tmat <- array(unlist(mypred), dim = mydim)
  colnames(res$Tmat) <- rownames(res$Tmat) <- names(object$prop)
  res$lags <- lags
  res$type <- "Theoretical"
  class(res) <- "transiogram"
  return(res)
}
