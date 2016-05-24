# Copyright (c) 2015 Santiago Barreda
# All rights reserved.

print.anova.rcr <-
function (x, ...){
  cat ("\nSignificance Tests for Groups of Coefficients\n\n")
  cat ("\nCall:\n\n")
  print (x$call)
  
  cat ("\n\n")
  printCoefmat (x$coefficients, has.Pvalue = TRUE)
}
