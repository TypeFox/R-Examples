library(bbmle)

x <- 0:10
y <- c(26, 17, 13, 12, 20, 5, 9, 8, 5, 4, 8)
d <- data.frame(x,y)
oldopts <- options(warn=-1,digits=3)  ## ignore warnings
m1 <- mle2(y~dpois(lambda=ymax/(1+x/xhalf)),
           start=list(ymax=1,xhalf=1),data=d)
m1
y2 <- c(26, 17, 10, 15, 20, 5, 9, 8, 5, 4, 8)
d2 <- data.frame(x,y=y2)

m2 <- update(m1,data=d2)
m2
m3 <- update(m1,.~dpois(lambda=c),start=list(c=5))
m3
options(oldopts)
