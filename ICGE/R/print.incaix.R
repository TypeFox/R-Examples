print.incaix <- function(x, ...){
# x must be a list with three components
   cat ("\n ---INCA index---\n")
   k <- length(x$well_class)
   cat("\n      number    well       %\n")
   cat("Group units   classified  well class. \n")
   for (i in 1:k){
       ni <- x$Ni_cluster[i]
       wi <- x$well_class[i]
       p <- wi/ni
       cat("  ", i,"      ", ni ,"      ", wi,"        ", format(p, digits=2), " \n", sep="")
   }
   cat("\n INCA index :", format(x$Total, digits=4), "\n")
   cat(" ---------------------------------------------\n")
}
