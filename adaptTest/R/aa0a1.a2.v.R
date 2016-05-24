`aa0a1.a2.v` <-
function (a,a0,a1, tol=.Machine$double.eps^.5)  ifelse(!(le(0,a1) && le(a1,a) && le(a,a0) && le(a0,1)),              NA, ifelse(eq(a0,a1),  1, ifelse(eq(a1,a), 0, ifelse(eq(a0,a), 1, ifelse(eq(a0,1)&&eq(a1,0), a, uniroot(function(x) a0a1a2.a.v(a0,a1,x) - a, lower=0, upper=1, tol=tol)$root)))))

