`a0a1a2.a.l` <-
function (a0,a1,a2) ifelse(!(le(0,a1) && le(a1,a0) && le(a0,1) && le(0,a2) && le(a2,1)), NA, ifelse(eq(a0,a1)||eq(a2,0), a1, ifelse(eq(a2,1), a0, ifelse(eq(a0,1)&&eq(a1,0), a2, a1 + int.decr(function (x) 1-pnorm(sqrt(2)*a2.c.l(a2)+qnorm(x)), a1, a0)))))

