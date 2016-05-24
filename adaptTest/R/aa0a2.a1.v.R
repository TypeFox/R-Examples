`aa0a2.a1.v` <-
function (a,a0,a2, tol=.Machine$double.eps^.5)  ifelse(!(le(0,a) && le(a,a0) && le(a0,1) && le(0,a2) && le(a2,1)),   NA, ifelse(eq(a,a0)||eq(a2,0), a, ifelse(eq(a,0)||eq(a2,1), NA, {ctmp <- a0a1a2.a.v(a0,0,a2); ifelse(gt(ctmp,a), NA, ifelse(eq(ctmp,a), 0, uniroot(function(x) a0a1a2.a.v(a0,x,a2) - a, lower=0, upper=a0, tol=tol)$root))})))

