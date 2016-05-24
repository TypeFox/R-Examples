print.glinternet = function(x, ...){

  sizes = t(sapply(x$activeSet[-1], function(x) sapply(x, function(y) ifelse(is.null(y), 0, nrow(y)))))
  sizes = rbind(0, sizes)
  output = data.frame(lambda=signif(x$lambda, 3), objValue=signif(x$objValue, 3), sizes)

  #print the call
  cat("Call: ")
  print(x$call)

  #print number of interactions found
  print(output)
}
