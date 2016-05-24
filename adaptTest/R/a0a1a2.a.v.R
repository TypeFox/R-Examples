`a0a1a2.a.v` <-
function (a0,a1,a2) ifelse(!(le(0,a1) && le(a1,a0) && le(a0,1) && le(0,a2) && le(a2,1)), NA, ifelse(eq(a0,a1)||eq(a2,0), a1, ifelse(eq(a2,1), a0, ifelse(eq(a0,1)&&eq(a1,0), a2, {ctmp <- a2.c.v(a2); a1 + int.decr(function(x) (1-x^ctmp)^(1/ctmp), a1, a0)}))))

