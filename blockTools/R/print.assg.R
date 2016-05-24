print.assg <- function(x, digits = max(2, getOption("digits") - 4),
                        ...){
  cat("\nAssignments:\n\n")
  if(length(x$assg) == 1){
    print.default(format(as.matrix(x$assg[[1]]), digits = digits),
                  print.gap = 2, quote = FALSE, ...)
  }else{
    for(i in 1:length(x$assg)){
      cat("Group:", names(x$assg)[i], "\n")
      print.default(format(as.matrix(x$assg[[i]]), digits = digits),
                    print.gap = 2, quote = FALSE, ...)
      cat("\n")
    }
  }
  invisible(x)
}
