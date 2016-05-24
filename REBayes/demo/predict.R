n <- 1000
m <- sample(c(0,4), n, prob = c(.7, .3), replace = TRUE) 
x <- rnorm(n,m)
f <- GLmix(x)
mse1 <- sum((predict(f,x,1)-m)^2)
mse2 <- sum((f$dy-m)^2)
print(c(mse1,mse2))
mae1 <- sum(abs(predict(f,x,1)-m))
mae2 <- sum(abs(f$dy-m))
print(c(mae1,mae2))


