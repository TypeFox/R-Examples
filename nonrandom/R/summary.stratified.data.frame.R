
summary.stratified.data.frame <- function(object,
                                          ...){

  df <-
    data.frame(cbind(object$intervals,
                     as.numeric(table(object$stratum.index)),
                     round(as.numeric(table(object$stratum.index))/dim(object$data)[1],3)))

  colnames(df) <- c("  Strata bounds","    n","    n (proportion)")

  sum.str <- list(str.var  = object$stratified.by,
                  str.info = df)

  class(sum.str) <- "summary.stratified.data.frame"

  sum.str

}


