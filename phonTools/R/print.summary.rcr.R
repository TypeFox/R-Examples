# Copyright (c) 2015 Santiago Barreda
# All rights reserved.


print.summary.rcr <-
function (x, ...){
  cat ("\nSignificance Tests for Individual Coefficients\n")
  cat ("\nCall:\n")
  print (x$call)
  
  cat ("\n")
  printCoefmat (x$coefficients, has.Pvalue = TRUE)
}
