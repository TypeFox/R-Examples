print.blockSum <- function(x, digits = max(2, getOption("digits") - 4),
                        ...){

  cat("Number of Groups =", length(x$blocks), "\n")
  if(length(x$blocks)==1){
    cat("Number of Blocks =", nrow(x$blocks[[1]]), "\n")
    cat("Mean distance per Block =", mean(x$blocks[[1]][,3]), "\n")
  }else{
    for(i in 1:length(x$blocks)){
      cat("Number of Blocks in Group", names(x$blocks)[i], "=",
          nrow(x$blocks[[i]]), "\n") 
      cat("Mean distance per Block, Group", names(x$blocks)[i], "=",
          mean(x$blocks[[i]][,3]), "\n")
    }
  }
  
#  cat("\nBlocks:\n\n")
#  if(length(x$blocks)==1){
#    print.default(format(as.matrix(x$blocks[[1]]), digits = digits),
#                  print.gap = 2, quote = FALSE, ...)
#  }else{
#    for(i in 1:length(x$blocks)){
#      cat("Group:", names(x$blocks)[i], "\n")
#      print.default(format(as.matrix(x$blocks[[i]]), digits = digits),
#                    print.gap = 2, quote = FALSE, ...)
#      cat("\n")
#    }
# }
  invisible(x)
}
