library(bbmle)
set.seed(1002)
lymax <- c(0,2)
lhalf <- 0
x <- runif(200)
g <- factor(rep(c("a","b"),each=100))
y <- rnbinom(200,mu=exp(lymax[g])/(1+x/exp(lhalf)),size=2)
d <- data.frame(x,g,y)

fit3 <- mle2(y~dnbinom(mu=exp(lymax)/(1+x/exp(lhalf)),size=exp(logk)),
    parameters=list(lymax~g),
    start=list(lymax=0,lhalf=0,logk=0),data=d)

plot(y~x,col=g)
## true curves
curve(exp(0)/(1+x/exp(0)),add=TRUE)
curve(exp(2)/(1+x/exp(0)),col=2,add=TRUE)
xvec = seq(0,1,length=100)
lines(xvec,predict(fit3,newdata=list(g=factor(rep("a",100),levels=c("a","b")),
                                x = xvec)),col=1,lty=2)
lines(xvec,predict(fit3,newdata=list(g=factor(rep("b",100),levels=c("a","b")),
                                x = xvec)),col=2,lty=2)

p1 = predict(fit3)
## manual prediction
p2A =
with(as.list(coef(fit3)),exp(`lymax.(Intercept)`)/(1+x[1:100]/exp(lhalf)))
p2B = with(as.list(coef(fit3)),exp(`lymax.(Intercept)`+lymax.gb)/(1+x[101:200]/exp(lhalf)))
p2 = c(p2A,p2B)
all(p1==p2)


