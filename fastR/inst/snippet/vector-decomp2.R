v1 <- c(1,0,0)
v2 <- c(1,1,1)
x  <- c( 2,3,5)
vdecomp(x,v1,v2)
h <- x - vdecomp(x,v1,v2)$x.proj; 
round(h,8)
round(dot(h,v1),8)
round(dot(h,v2),8)
