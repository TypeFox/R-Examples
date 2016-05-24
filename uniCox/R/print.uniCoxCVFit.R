print.uniCoxCVFit=function(x,...){
  cat("Call:\n")
  dput(x$call)
  mat <- rbind(number.of.genes = format(round(x$ncallcvm, 3)), ave.drop.in.deviance =round(x$devcvm,3), se=format(round(x$se.devcvm, 3))) 
  dimnames(mat) <- list(dimnames(mat)[[1]], paste(1:ncol(mat)))
  print(t(mat), quote = FALSE)
  invisible()
}

