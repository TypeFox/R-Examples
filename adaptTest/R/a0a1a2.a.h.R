`a0a1a2.a.h` <-
function (a0,a1,a2) ifelse(!(le(0,a1) && le(a1,a0) && le(a0,1) && le(0,a2) && le(a2,1)), NA, a1 + a2*(a0-a1))

