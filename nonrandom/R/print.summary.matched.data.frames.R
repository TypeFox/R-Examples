
print.summary.matched.data.frames <- function(x,
                                              ...){

  cat("\n Matched by: ", x$match.var, "\n")
  cat("\n Matching parameter:\n")
  print(x$parameter)
  cat("\n Matching information:\n")
  print(x$info)
  cat("\n Matching data:\n")
  print(x$number)

}


