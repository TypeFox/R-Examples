print.summary.progenyClust <-
function(x,...){
  cat("Call:\n")
  print(x$call)
  cat('\n')
  cat("Optimal Number of Clusters:\n")
  if(!is.na(x$n.gap)){
    cat(paste0("gap criterion - ",x$n.gap))
  }
  if(!is.na(x$n.score)){
    cat(paste0("\nscore criterion - "),x$n.score)
  }
  cat(' \n\n')
}
