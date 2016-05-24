`a2.c.b` <-
function (a2)       ifelse(!(le(0,a2) && le(a2,1)),                                      NA, exp(-qchisq(1-a2,4)/2))

