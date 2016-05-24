ciTest <- function(x,set=NULL,...){
  UseMethod("ciTest")
}

ciTest.table <- function(x, set=NULL, ...){
  ciTest_table(x, set, ...)
}

ciTest.list <- function(x, set=NULL, ...){
  ciTest_mvn(x, set, ...)
}

ciTest.data.frame <- function(x, set=NULL, ...){
  ciTest_df(x, set, ...)
}


print.citest <- function(x,...){
  if (length(x$varNames) > 2){
    cat("Testing", x$varNames[1], "_|_", x$varNames[2], "|",x$varNames[-(1:2)],"\n")
  } else {
    cat("Testing", x$varNames[1], "_|_", x$varNames[2], "\n")
  }
  cat(sprintf("Statistic (%s): %8.3f df: %s p-value: %6.4f method: %s\n",
              x$statname, x$statistic, x$df, x$p.value, x$method))
}

summary.citest <- function(object,...){
  str(object,max.level=1)
  return(invisible(object))
}
