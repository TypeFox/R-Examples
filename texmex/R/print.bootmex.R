print.bootmex <- function(x , ... ){
  print(x$call )
  cat(paste( "\n", x$R, " bootstrap samples created.\n\n" , sep = "" ) )

  # Sometimes bootstrap samples fail to get convergence of the optimizer.
  # Check to see how many effective bootstrap samples were performed.
  object <- x$boot
  wh <- unlist(lapply(object, function(z) dim(z$Z)[[1]] == dim(na.omit(z$Z))[[1]]))
  object <- object[wh]
  eff <- sum(wh)

  cat( paste("Dependence structure estimation successful for", eff, "effective samples.\n" ) )
#  cat( paste("Dependence structure bootstrap mean parameter estimates:\n" ))
  cat("\n")
  invisible()
}