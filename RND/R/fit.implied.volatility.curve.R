fit.implied.volatility.curve <-
function(x,k)
{

  obj = lm(x ~ k + I(k^2) )
  a0 = obj$coefficients[[1]]
  a1 = obj$coefficients[[2]]
  a2 = obj$coefficients[[3]]
  
  summary.obj = summary(obj)
  out = list(a0 = a0, a1 = a1, a2 = a2, summary.obj = summary.obj)
  out
}
