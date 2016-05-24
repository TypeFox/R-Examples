girls <- subset(juul2, age<20 & age>5 & sex==2)
boys <- subset(juul2, age<20 & age>5 & sex==1)
young <- subset(juul2, age<20 & age>5)
stval <- c(alpha=exp(5.3),beta=exp(0.42),gamma=0.15)
fit.boys <- nls(height~alpha*exp(-beta*exp(-gamma*age)),
       start=stval, data=boys)
fit.girls <- nls(height~alpha*exp(-beta*exp(-gamma*age)),
       start=stval, data=girls)
fit.young <- nls(height~alpha*exp(-beta*exp(-gamma*age)),
       start=stval, data=young)
ms.pooled <- (deviance(fit.boys) + deviance(fit.girls))/(499+625)
ms.diff <- (deviance(fit.young) - 
            deviance(fit.boys) - deviance(fit.girls))/3
ms.diff/ms.pooled   
fit.young2 <- nls(height~(alpha+da*(sex==1))*
                   exp(-(beta+db*(sex==1))*
                   exp(-(gamma+dg*(sex==1))*age)),
           start=c(alpha=exp(5.3),beta=exp(0.42),gamma=0.15,
           da=0, db=0, dg=0), data=young)
summary(fit.young2)
anova(fit.young, fit.young2)
e1 <- subset(philion, experiment==1)
fit <- nls(sqrt(response) ~ sqrt(ymax / (1 + (dose/ec50)^exp(la))), 
        start=list(ymax=28, ec50=.3, la=0), data=e1, 
        lower=c(.1,.0001,-Inf), algorithm="port")
summary(fit)
confint(fit)
p <- profile(fit, alphamax=.2)
par(mfrow=c(3,1))
plot(p)
confint(p)
e1 <- subset(philion, experiment==1)
fit1 <- nls(sqrt(response) ~ sqrt(ymax / (1 + dose/b)^exp(la)), 
        start=list(ymax=28, b=.3, la=0), data=e1, 
        lower=c(.1,.0001,-Inf), algorithm="port")
summary(fit1)
fit2 <- nls(sqrt(response) ~ sqrt(ymax / (1 +
        dose/(d50/(2^(1/exp(la))-1)))^exp(la)),
        start=list(ymax=28, d50=.3, la=0), data=e1,
        lower=c(.1,.0001,-Inf), algorithm="port")
summary(fit2)
dd <- seq(0,1,,200)
yy <- predict(fit, newdata=data.frame(dose=dd))
y1 <- predict(fit2, newdata=data.frame(dose=dd))
matplot(dd,cbind(yy,y1)^2, type="l")  
rm(list=ls())
while(search()[2] != "package:ISwR") detach()
