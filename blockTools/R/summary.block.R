summary.block.R <- function(object, ...){
#  if(length(object$blocks)==1){
#    cat("Number of Groups = 1\n")
#    cat("Number of Blocks =", nrow(x$blocks[[1]]), "\n")
#    cat("Mean distance per Block =", mean(x$blocks[[1]][,3]), "\n")
#  object$blocks <- c(1, nrow(object$blocks[[1]]), mean(object$blocks[[1]][,3]))
#  colnames(object$blocks) <- c("Estimate")
#  rownames(object$blocks) <- c("Number of Groups", "bl", "dis")  
  class(object) <- "summary.block"
  object
}


 #   for(i in 1:length(x$blocks))




  #object$short <- short
  #class(object) <- "eiMDsum"
  #object
