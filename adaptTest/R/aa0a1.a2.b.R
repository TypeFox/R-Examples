`aa0a1.a2.b` <-
function (a,a0,a1, tol=.Machine$double.eps^.5)  ifelse(!(le(0,a1) && le(a1,a) && le(a,a0) && le(a0,1)),              NA, ifelse(eq(a0,a1),  1, ifelse(eq(a1,a), 0, ifelse(eq(a0,a), 1, ifelse(eq(a0,1), ifelse(le(a1,a2.c.b(a)), a, c.a2.b((a1-a)/log(a1))), {ctmp <- uniroot(function(x)x*(1+log(a0)-log(x))-a, lower=a2.c.b(a), upper=a0, tol=tol)$root; ifelse(le(a1,ctmp), c.a2.b(ctmp), c.a2.b((a-a1)/(log(a0)-log(a1))))})))))

