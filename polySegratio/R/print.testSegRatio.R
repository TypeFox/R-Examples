`print.testSegRatio` <-
function(x, ..., last=10)
{
  cat("Segregation ratio test:\n")
  cat("Method:",x$method,"at the",x$alpha,"level\n")
  cat("Expected segregation ratios:\n")
  print(unlist(x$E.segRatio), ...)
  
  cat("Proportion of markers classified at",x$alpha,"level:",
      sum(rowSums(x$allocated)==1)/dim(x$allocated)[1],"\n")
  cat("Classified:", sum(rowSums(x$allocated)==1),
      ", not classified:", 
      sum(rowSums(x$allocated)!=1),"\n")

  ## classified twice (includes all classified as 4 dosage)
  markers.dc <- x$allocated[rowSums(x$allocated)>1,]
  no.dc <- sum(rowSums(x$allocated)>1)
  if (length(markers.dc)<1) {
    cat("No markers doubly classified\n")
  } else {
    cat("Markers doubly classified:",no.dc,"\n")
    if (no.dc==1) cat(rownames(x$allocated)[rowSums(x$allocated)>1],"\n")
    print(markers.dc)
  }

  cat("Number in each marker dosage class (classified once):")
  print(table(factor(x$dosage,
                     labels=names(x$E.segRatio$ratio),
                     levels=1:length(x$E.segRatio$ratio))), ...)
  cat("Dosage of first",last,"markers (where dosage unique):\n")
  print(x$dosage[1:last], ...)
  cat("Call: ")
  print(x$call)

}

