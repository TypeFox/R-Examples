get.risk.t0 <-
function(coef, marker, linkinvfun) linkinvfun(coef[1] + coef[3]*marker)
