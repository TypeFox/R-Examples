`aa0a2.a1.h` <-
function (a,a0,a2)  ifelse(!(le(0,a) && le(a,a0) && le(a0,1) && le(0,a2) && le(a2,1)),   NA, ifelse(eq(a,a0)||eq(a2,0), a, ifelse(eq(a,0)||eq(a2,1), NA, ifelse(gt(a2*a0,a), NA, (a-a2*a0)/(1-a2)))))

