# Copyright (c) 2015 Santiago Barreda
# All rights reserved.


summary.rcr <-
function (object, ...){
  if (sum(is.na(object$coefficients)) > 0) return (cat ("Error: Null coefficient values. Check individual coefficients and model fits.\n\n"))
  n = nrow (object$coefficients)
  estimate = colMeans (object$coefficients)
  se = sqrt(diag(var(object$coefficients))) / sqrt(n)
  t.value = estimate / se
  p.value = pt (-abs(t.value), n-1)*2 

  coefficients =  cbind (estimate = estimate, se = se, t.value = t.value, df = n-1, p.value = p.value)
  output = list (call = object$call, coefficients = coefficients)
  class(output) = "summary.rcr"
  output
}
