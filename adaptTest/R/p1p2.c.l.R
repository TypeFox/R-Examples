`p1p2.c.l` <-
function (p1,p2)    ifelse(!(le(0,p1) && le(p1,1) && le(0,p2) && le(p2,1)),              NA, ifelse((eq(p1,0)&&gt(p2,0))||(eq(p1,1)&&lt(p2,1)), NA, (qnorm(1-p1)+qnorm(1-p2))/sqrt(2)))

