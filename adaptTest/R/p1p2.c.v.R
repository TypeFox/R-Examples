`p1p2.c.v` <-
function (p1,p2, tol=.Machine$double.eps^.5)    ifelse(!(le(0,p1) && le(p1,1) && le(0,p2) && le(p2,1)),              NA, ifelse((eq(p1,0)&&gt(p2,0))||(eq(p1,1)&&lt(p2,1)), NA, ifelse(lt(p1,1)&&eq(p2,0), 0, ifelse(gt(p1,0)&&eq(p2,1), Inf, uniroot(function(x) p1^x+p2^x - 1, lower=0, upper=10^20, tol=tol)$root))))

