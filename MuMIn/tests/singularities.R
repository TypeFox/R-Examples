library(MuMIn)
options(na.action = "na.fail")

set.seed(1)
zz <- data.frame(x=runif(15), f1=gl(3,5), f2=factor(rep(1:2,c(10,5))))
zz$y <- 100*zz$x + as.numeric(zz$f1)*10 * as.numeric(zz$f2)

nafit <- lm(y~f1*f2*x, zz)

summary(nafit)
coef(nafit)

gm <- get.models(dredge(nafit), subset = NA)
ma <- model.avg(gm)

summary(ma)
coef(ma, T)
confint(ma)

predict(ma)

#Sys.sleep(5)
