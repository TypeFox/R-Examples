suppressMessages(library(cobs))

options(digits = 6)

x <- c(1:3, 5,6,9,12)
y <- c(-1:1,0,1,-2,0) + 8*x

all(abs(residuals(cobs(x,y))) < 1e-10)# no constraints -> perfect fit
summary(c1 <- cobs(x,y, constraint = "convex"))#  99.88 %
summary(c2 <- cobs(x,y, constraint = "concave"))# 99.84 %

sum(resid(c1)^2)# 7
sum(resid(c2)^2)# 9.715

plot(x,y)
lines(predict(c1), col=2)
lines(predict(c2), col=3)

##using browser() inside
predict(c1, interval = "both")
predict(c2, interval = "both")

cat('Time elapsed: ', proc.time(),'\n') # for ``statistical reasons''
