`aa1a2.a0.b` <-
function (a,a1,a2, tol=.Machine$double.eps^.5)  ifelse(!(le(0,a1) && le(a1,a) && le(a,1) && le(0,a2) && le(a2,1)),   NA, ifelse(eq(a,a1)||eq(a2,1), a, ifelse(eq(a,1)||eq(a2,0), NA, ifelse(ge(a2.c.b(a2),a), a, {ctmp <- a0a1a2.a.b(1,a1,a2); ifelse(lt(ctmp,a), NA, ifelse(eq(ctmp,a), 1, uniroot(function(x) a0a1a2.a.b(x,a1,a2) - a, lower=a1, upper=1, tol=tol)$root))}))))

