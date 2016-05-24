`a0a1a2.a.b` <-
function (a0,a1,a2) ifelse(!(le(0,a1) && le(a1,a0) && le(a0,1) && le(0,a2) && le(a2,1)), NA, ifelse(eq(a0,a1)||eq(a2,0), a1, ifelse(eq(a2,1), a0, {ctmp <- a2.c.b(a2); ifelse(ge(ctmp,a0), a0, ifelse(ge(ctmp,a1), ifelse(eq(a0,1), a2, ctmp*(1+log(a0)-log(ctmp))), a1+ctmp*(log(a0)-log(a1))))})))

