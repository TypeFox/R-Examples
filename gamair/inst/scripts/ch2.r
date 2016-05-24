# R code for chapter 2 of Wood (2006) "GAMs: An Introduction with R"

## 2.2.2 Geometry and IRLS convergence

x <- c(.6,1.5);y <- c(.02,.9)
ms <- exp(-x*4)   # set initial values at lower left
glm(y~I(-x)-1,family=gaussian(link=log),mustart=ms)
ms <- exp(-x*0.1)  # set initial values at upper right
glm(y~I(-x)-1,family=gaussian(link=log),mustart=ms)

## 2.3.1 Binomial models and heart disease

par(mfrow=c(1,1))
heart <- data.frame(ck = 0:11*40+20,
ha=c(2,13,30,30,21,19,18,13,19,15,7,8),
ok=c(88,26,8,5,0,1,1,1,1,0,0,0))
p<-heart$ha/(heart$ha+heart$ok)
plot(heart$ck,p,xlab="Creatinine kinase level",
     ylab="Proportion Heart Attack")

mod.0 <- glm(cbind(ha,ok)~ck,family=binomial(link=logit),data=heart)
mod.0 <- glm(cbind(ha,ok)~ck,family=binomial,data=heart)

mod.0

(271.7-36.93)/271.7

1-pchisq(36.93,10)

op <- par(mfrow=c(2,2))
plot(mod.0)

par(op)
plot(heart$ck,p,xlab="Creatinine kinase level",ylab="Proportion Heart Attack")
lines(heart$ck,fitted(mod.0))

mod.2 <- glm(cbind(ha,ok)~ck+I(ck^2)+I(ck^3),family=binomial,data=heart)
mod.2

par(mfrow=c(2,2))
plot(mod.2)

par(mfrow=c(1,1))
plot(heart$ck,p,xlab="Creatinine kinase level",ylab="Proportion Heart Attack")
lines(heart$ck,fitted(mod.2))

anova(mod.0,mod.2,test="Chisq")

## 2.3.2 A Poisson regression epidemic model

y <- c(12,14,33,50,67,74,123,141,165,204,253,246,240)
t <- 1:13
plot(t+1980,y,xlab="Year",ylab="New AIDS cases",ylim=c(0,280))

m0 <- glm(y~t,poisson)
m0
par(mfrow=c(2,2))
plot(m0)

m1 <- glm(y~t+I(t^2),poisson)
plot(m1)
summary(m1)

anova(m0,m1,test="Chisq")

beta.1 <- summary(m1)$coefficients[2,]
ci <- c(beta.1[1]-1.96*beta.1[2],beta.1[1]+1.96*beta.1[2])
ci # print 95% CI for beta_1


par(mfrow=c(1,1))
new.t<-seq(1,13,length=100)
fv <- predict(m1,data.frame(t=new.t),se=TRUE)
plot(t+1980,y,xlab="Year",ylab="New AIDS cases",ylim=c(0,280))
lines(new.t+1980,exp(fv$fit))
lines(new.t+1980,exp(fv$fit+2*fv$se.fit),lty=2)
lines(new.t+1980,exp(fv$fit-2*fv$se.fit),lty=2)

## 2.3.3 Log-linear models for catagorical data

al <- data.frame(y=c(435,147,375,134),gender=as.factor(c("F","F","M","M")),
                 faith=as.factor(c(1,0,1,0)))
al

mod.0 <- glm(y~gender+faith,data=al,family=poisson)
model.matrix(mod.0)

mod.0
fitted(mod.0)

mod.1<-glm(y~gender*faith,data=al,family=poisson)
model.matrix(mod.1)
mod.1

anova(mod.0,mod.1,test="Chisq")

## 2.3.4 Sole eggs in the Bristol channel

library(gamair)
data(sole)
sole$off <- log(sole$a.1-sole$a.0)# model offset term
sole$a<-(sole$a.1+sole$a.0)/2     # mean stage age
solr<-sole                        # make copy for rescaling
solr$t<-solr$t-mean(sole$t)
solr$t<-solr$t/var(sole$t)^0.5
solr$la<-solr$la-mean(sole$la)
solr$lo<-solr$lo-mean(sole$lo)
b <- glm(eggs ~ offset(off)+lo+la+t+I(lo*la)+I(lo^2)+I(la^2)
               +I(t^2)+I(lo*t)+I(la*t)+I(lo^3)+I(la^3)+I(t^3)
               +I(lo*la*t)+I(lo^2*la)+I(lo*la^2)+I(lo^2*t)
               +I(la^2*t)+I(la*t^2)+I(lo*t^2)+ a +I(a*t)+I(t^2*a),
         family=quasi(link=log,variance="mu"),data=solr)
summary(b)

b1<-update(b,~.-I(lo*t))
summary(b1) # better plot

b2<-update(b1,~.-I(lo*la*t))
summary(b2) # worse again

b3<-update(b2,~.-I(lo*t^2))
summary(b3) # better again

b4<-update(b3,~.-I(lo^2*t))
summary(b4) # better again

anova(b,b4,test="F")

par(mfrow=c(1,2)) # split graph window into 2 panels
plot(fitted(b4)^0.5,solr$eggs^0.5) # fitted vs. data plot
plot(fitted(b4)^0.5,residuals(b4)) # resids vs. sqrt(fitted)
