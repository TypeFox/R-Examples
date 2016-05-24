precision <- function(X, Hz){
  end <- length(X) / (Hz / 10)
  window <- rep(1:end, each = (Hz / 10))
  Smooth <- as.vector(unlist(by(X[1:length(window)], window, MFW)))
  
  return(mean(abs(X[1:length(window)] - Smooth), na.rm = T))
}
