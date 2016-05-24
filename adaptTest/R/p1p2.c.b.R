`p1p2.c.b` <-
function (p1,p2)    ifelse(!(le(0,p1) && le(p1,1) && le(0,p2) && le(p2,1)),              NA, ifelse(eq(p1,0)&&gt(p2,0), NA, p1*p2))

