print.summary.BANOVA.Multinomial <-
function(x, ...){
  cat('Call:\n')
  print(x$call)
  
  cat('\nConvergence diagnostics:\n')
  print(x$conv)
  
  cat('\nTable of sum of squares and effect sizes (Bayesian ANOVA/ANCOVA):\n')
  print(x$anova.table)
  
  cat('\nTable of p-values (Multidimensional): \n')
  print(x$pvalue.table)
  
  cat('\nTable of coefficients: \n')
  printCoefmat(x$coef.table)
  
}
