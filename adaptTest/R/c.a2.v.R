`c.a2.v` <-
function (c)        ifelse(!(le(0,c)),                                                   NA, ifelse(eq(c,0), 0, ifelse(eq(c,Inf), 1, exp(lgamma(1/c)-lgamma(1/c+.5)-2*log(2)/c)*sqrt(pi)/c)))

