`print.runJagsWrapper` <-
function(x, ...){
  
  if (class(x) != "runJagsWrapper")
    stop("'x' must be of class 'runJagsWrapper'")

  print(x$run.jags)
  print(x$summary)
  print(x$diagnostics)
  print(x$doses)
  cat("\nDIC:",x$DIC,"\n")
}

