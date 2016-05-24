v1 <- c(1,0,0)
v2 <- c(1,1,1)
v3 <- c(1,2,3)
x  <- c(2,7,3)
vdecomp(x,v1,v2)
a1 <- vdecomp(x,v1,v2)$coef[1]
b1 <- vdecomp(x,v1,v2)$coef[2]
x1 <- a1 * v1 + b1 * v2; x1
# decompose x into x1 and v3 
vdecomp(x,x1,v3)
a2 <- vdecomp(x,x1,v3)$coef[1]
b2 <- vdecomp(x,x1,v3)$coef[2]
# this should equal x
a2 * (a1 * v1 + b1* v2) + b2 * v3 
# the three coefficients
c( a2 * a1, a2 * b1, b2)
