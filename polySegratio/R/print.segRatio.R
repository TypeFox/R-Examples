`print.segRatio` <-
function(x, digits=3, ..., index=c(1:min(10,length(x$r))) ) {
  
  cat("Summary statistics for segregation ratios:\n")
  print(summary(x$seg.ratio),...)

  cat("Observed numbers and segregation proportions for\n",
      length(index),"of the markers for",x$n.individuals,
      "individuals:\n")
  miss <- x$n.individuals*length(x$n) - sum(x$n)
  if( miss>0 ) {
    cat("Percentage of missing markers:",100*miss/sum(x$n),"\n")
  }
  ## print(cbind(r=x$r[index], n=x$n[index],
  ##       segRatio=x$x[index]))
  print(x$seg.ratio[index],digits=digits, ...)

}

