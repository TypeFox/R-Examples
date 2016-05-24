# Copyright (c) 2015 Santiago Barreda
# All rights reserved.


anova.rcr = function (object, ...){
  if (sum(is.na(object$coefficients)) > 0) return (cat ("Error: Null coefficient values. Check individual coefficients and model fits.\n\n"))
  
  nfactors = length (levels (object$factors))
  factors = levels (object$factors)
  nparts = length (levels(object$participants))

  f.value = NULL; df1 = NULL; df2 = NULL;

  for (i in 1:nfactors){
    temp = object$coefficients[,object$factors == factors[i]]
    if (is.null(ncol(temp))){
      f.value = c(f.value, t.test(temp)$statistic^2)
      df1 = c(df1, 1)
      df2 = c(df2, nparts - 1)
    }  
    else if (!is.null(ncol(temp))){
      f.value = c(f.value, hotelling.test (temp)$f.value)
      df1 = c(df1, ncol (temp))
      df2 = c(df2, nparts - ncol (temp))
    }
  }
  p.value = 1 - pf (f.value, df1, df2)
  coefficients = data.frame (df1 = df1, df2 = df2, f.value = f.value, p.value = p.value)
  rownames (coefficients) = object$factor.names

  output = list (call = object$call, coefficients = coefficients)

  class(output) = "anova.rcr"
  output
}


