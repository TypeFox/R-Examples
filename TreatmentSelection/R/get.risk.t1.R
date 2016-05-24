get.risk.t1 <-
function(coef, marker, linkinvfun) linkinvfun(coef[1] + coef[2] + coef[3]*marker + coef[4]*marker)
