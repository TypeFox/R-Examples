`c.a2.b` <-
function (c)        ifelse(!(le(0,c) && le(c,1)),                                        NA, 1-pchisq(-2*log(c), 4))

