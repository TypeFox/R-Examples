print.pendensity <- function(x,...) {
  obj <- x
  K.help <- (obj$splines$K-1)/2
  K.seq <- -K.help:K.help
  cat("\nCall:",deparse(obj$call),"\n\n", sep="")
  cat("\nWeights:\n",sep="")
  mat <- matrix(obj$results$ck,obj$splines$N,obj$splines$K)
  colnames(mat) <- paste("K=",K.seq,sep="")
  print(mat)
  cat("\nlambda0:",obj$results$lambda0,"\n",sep="")
  cat("\nAIC:",obj$results$AIC$my.AIC,"\n",sep="")
}