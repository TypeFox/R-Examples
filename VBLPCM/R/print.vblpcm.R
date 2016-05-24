print.vblpcm<-function(x, ...)
  {
  cat("Variational Bayes Approximation to Latent Position Cluster Model\n")
  cat("Number of nodes:\t",x$N,"\n")
  cat("Number of edges:\t",x$NE,"\n")
  cat("Number of groups:\t",x$G,"\n")
  cat("Latent Space dimension:\t",x$d,"\n")
  cat("Model:\t", x$model, "\n")
  }
