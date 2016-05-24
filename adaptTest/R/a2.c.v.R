`a2.c.v` <-
function (a2, tol=.Machine$double.eps^.5)       ifelse(!(le(0,a2) && le(a2,1)),                                      NA, ifelse(eq(a2,0), 0, ifelse(eq(a2,1), Inf, uniroot(function(x) c.a2.v(x) - a2, lower=0, upper=10^10, tol=tol)$root)))

