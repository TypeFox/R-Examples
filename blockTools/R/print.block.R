print.block <- function(x, digits = max(2, getOption("digits") - 4),
                        ...){
  cat("\nBlocks:\n\n")
  if(length(x$blocks) == 1){
    print.default(format(as.matrix(x$blocks[[1]]), digits = digits),
                  print.gap = 2, quote = FALSE, ...)
  }else{
    for(i in 1:length(x$blocks)){
      cat("Group:", names(x$blocks)[i], "\n")
      print.default(format(as.matrix(x$blocks[[i]]), digits = digits),
                    print.gap = 2, quote = FALSE, ...)
      cat("\n")
    }
  }
  invisible(x)
}
