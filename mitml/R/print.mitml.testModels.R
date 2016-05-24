print.mitml.testModels <- function(x,...){
# print method for MI estimates

  cl <- x$call
  test <- x$test
  mth <- x$method
  use <- x$use
  reml <- x$reml
  m <- x$m
  adj <- x$adj.df
  dfc <- x$df.com

  # header
  cat("\nCall:\n", paste(deparse(cl)), sep="\n")

  cat("\nModel comparison calculated from",m,"imputed data sets.")
  cat("\nCombination method:",mth,
    if(mth=="D2"){paste("(",use,")",sep="")},"\n")

  # check for large values
  fmt <- c("%.3f","%.0f","%.3f","%.3f","%.3f")
  fmt[test>=10^5] <- "%.3e"
  out <- sprintf(fmt,test)

  # print table
  cat("\n")
  w <- max(sapply(c(out,colnames(test)),nchar))
  cat("  ",format(colnames(test),justify="right",width=w),"\n")
  cat("  ",format(out,justify="right",width=w),"\n")

  if(mth=="D1"){
  cat(if(adj){c("\nHypothesis test adjusted for small samples with",
              paste("df=[",paste(dfc,collapse=","),"]\ncomplete-data degrees of freedom.",sep=""))
      }else{"\nUnadjusted hypothesis test as appropriate in larger samples."},"\n")
  }
  if(reml){
  cat("\nModels originally fit with REML were automatically refit using ML.\n")
  }

  cat("\n")

  invisible()
}

summary.mitml.testModels <- function(object,...){
# summary method for objects of class mitml.testModels

  print.mitml.testModels(object,...)

}
