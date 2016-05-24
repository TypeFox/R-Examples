`a2.c.l` <-
function (a2)       ifelse(!(le(0,a2) && le(a2,1)),                                      NA, qnorm(1-a2))

