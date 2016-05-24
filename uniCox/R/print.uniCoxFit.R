print.uniCoxFit=function(x,...){
  cat("Call:\n")
  dput(x$call)
  mat <- rbind(lambda = format(round(x$lamlist, 3)), number.of.genes =x$nfeatures) 
  dimnames(mat) <- list(dimnames(mat)[[1]], paste(1:ncol(mat)))
  print(t(mat), quote = FALSE)
  invisible()
}

