print.diagnose <- function(x, digits = max(2, getOption("digits") - 4),
                        ...){
  cat("\nUnits differing by at least ", x$suspect.range[1], " and no
more than ", x$suspect.range[2], " on ", x$suspect.var, ":\n\n", sep="")
  if(length(x$diagnose) == 1){
    print.default(format(as.matrix(x$diagnose[[1]]), digits = digits),
                  print.gap = 2, quote = FALSE, ...)
  }else{
    for(i in 1:length(x$diagnose)){
      cat("Group:", names(x$diagnose)[i], "\n")
      print.default(format(as.matrix(x$diagnose[[i]]), digits = digits),
                    print.gap = 2, quote = FALSE, ...)
      cat("\n")
    }
  }
  invisible(x)
}
