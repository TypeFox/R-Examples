
print.summary.stratified.pscore <- function(x,
                                            ...){

  cat("\n Stratified by: ", x$str.var, "\n")
  cat("\n Strata information:\n")
  print(x$str.info)

}


