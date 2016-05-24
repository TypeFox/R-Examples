`aa1a2.a0.h` <-
function (a,a1,a2)  ifelse(!(le(0,a1) && le(a1,a) && le(a,1) && le(0,a2) && le(a2,1)),   NA, ifelse(eq(a,a1)||eq(a2,1), a, ifelse(eq(a,1)||eq(a2,0), NA, ifelse(lt(a1+a2*(1-a1),a), NA, a1+(a-a1)/a2))))

