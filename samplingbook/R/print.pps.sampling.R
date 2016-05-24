print.pps.sampling <-
function(x,...){
  cat("\npps.sampling object: Sample with probabilities proportional to size\n")
  if(x$call$method=="sampford"){  cat("Method of Sampford:\n\n") } 
  if(x$call$method=="tille"){  cat("Method of Tille:\n\n") }
  if(x$call$method=="midzuno"){  cat("Method of Midzuno:\n\n") }    
  if(x$call$method=="madow"){  cat("Method of Madow:\n\n") } 
  cat("PPS sample: \n",sep="")
  print(x$sample)
  cat("\nSample probabilities: \n",sep="")
  print(x$PI)
  if (! is.null(x$PI.full)){
    cat("\nPairwise inclusion probabilities: \n",sep="")
    print(x$PI.full)
  } 
  invisible(x) 
}
