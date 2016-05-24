print.threshold <- function(x,...)
{
    cat("\n")
    for (i in 1:length(x$threshtable)) {
      cat("Design Matrix Block ",i,":\n",sep="")
      print(round(x$threshtable[[i]],5))
      cat("\n")
    }
    invisible(x$threshtable)
}