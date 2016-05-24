require(bbmle)
x <- 0:10
y <- c(26, 17, 13, 12, 20, 5, 9, 8, 5, 4, 8)
d <- data.frame(x,y)
fit <- mle2(y~dpois(lambda=ymax/(1+x/xhalf)), start=list(ymax=25,xhalf=3),data=d)
fit2 <- mle2(y~dpois(lambda=(x+1)*slope), start=list(slope=1),data=d)
BIC(fit)
BIC(fit,fit2)
