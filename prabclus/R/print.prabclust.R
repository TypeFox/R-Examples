"print.prabclust" <-
function(x, bic=FALSE, ...){
  cat("* Clustered presence-absence matrix * \n\n")
  cat("Clustered: ", x$mdsdim, "-dim. MDS result from method ",x$mdsmethod,"\n\n")
#  cat("mclust decided for ",ncol(x$clustsummary[[3]]$mu),
#      " clusters plus noise.\n")
  cat("Noise-detector NNclean has been used with k=", x$nnk, "\n")
  cat("NNclean is explained in S. Byers and A. E. Raftery, JASA 95 (1998), 781-794\n")
  cat("A Normal mixture model with noise component (mclust) has been used.\n")
  if (bic){
    cat("Cluster summary:\n")
    print(x$clustsummary)
    cat("\n")
    cat("BIC value summary:\n")
    print(x$bicsummary)
    cat("\n")
  }
  else{
    cat("Mixture component memberships:\n")
    print(x$clustering)
    cat("\n")
  }
  cat("Clustering (N denotes noise or one-point components):\n")
  print(x$symbols)
}
