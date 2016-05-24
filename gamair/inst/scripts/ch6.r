# R code for chapter 6 of Wood (2006) "GAMs: An Introduction with R"

library(gamair)
data(stomata)

m1 <- lm(area ~ CO2 + tree,stomata)
m0 <- lm(area ~ CO2,stomata)
anova(m0,m1)

m2 <- lm(area ~ tree,stomata)
anova(m2,m1)

st <- aggregate(data.matrix(stomata),
                by=list(tree=stomata$tree),mean)

st$CO2 <- as.factor(st$CO2);st

m3 <- lm(area~CO2,st)
anova(m3)

summary(m3)$sigma^2 - summary(m1)$sigma^2/4

## 6.1.3 A single random factor

library(nlme)  # load nlme `library', which contains data
data(Rail)     # load data
Rail
m1 <- lm(travel ~ Rail,Rail)
anova(m1)


# average over Rail effect...
rt <- aggregate(data.matrix(Rail),by=list(Rail$Rail),mean)
rt

m0 <- lm(travel ~ 1,rt)   # fit model to aggregated data
sigb <- (summary(m0)$sigma^2-summary(m1)$sigma^2/3)^0.5
# sigb^2 is variance component for rail
sig <- summary(m1)$sigma # sig^2 is resid. var. component
sigb
sig

summary(m0)

## 6.1.4 A model with two factors

library(nlme)
data(Machines)
names(Machines)
attach(Machines)  # make data available without `Machines$'
interaction.plot(Machine,Worker,score)

m1 <- lm(score ~ Worker*Machine,Machines)
m0 <- lm(score ~ Worker + Machine,Machines)
anova(m0,m1)

summary(m1)$sigma^2

Mach <- aggregate(data.matrix(Machines),by=
        list(Machines$Worker,Machines$Machine),mean)
Mach$Worker <- as.factor(Mach$Worker)
Mach$Machine <- as.factor(Mach$Machine)

m0 <- lm(score ~ Worker + Machine,Mach)
anova(m0)

summary(m0)$sigma^2 - summary(m1)$sigma^2/3

M <- aggregate(data.matrix(Mach),by=list(Mach$Worker),mean)
m00 <- lm(score ~ 1,M)
summary(m00)$sigma^2 - (summary(m0)$sigma^2)/3

## 6.2.2 Directly maximizing a mixed model likelihood in R

ll <- function(gamma,X,Z,y)
{ # untransform parameters
  sigma.b <- exp(gamma[1])
  sigma <- exp(gamma[2])
  n <- nrow(Z)

  # evaluate covariance matrix for y
  V <- Z%*%t(Z)*sigma.b^2 + diag(n)*sigma^2
  L <- chol(V) # L'L=V

  # transform dependent linear model to indep.
  y <- backsolve(L,y,transpose=TRUE)
  X <- backsolve(L,X,transpose=TRUE)
  b <- coef(lm(y~X-1)) # estimate fixed effects

  # evaluate log likelihood
  logLik <- -n/2*log(2*pi) -  sum (log(diag(L))) -
             sum((y-X%*%b)^2)/2
  attr(logLik,"fixed") <- b # allow retrieval of beta
  logLik
}

options(contrasts=c("contr.treatment","contr.treatment"))
Z <- model.matrix(~Rail$Rail-1)
X <- matrix(1,18,1)

rail.mod <- optim(c(0,0),ll,control=list(fnscale=-1),
                    X=X,Z=Z,y=Rail$travel)

exp(rail.mod$par) # variance components
attr(ll(rail.mod$par,X,Z,Rail$travel),"fixed")

## 6.3 Linear mixed models in R

library(nlme)
lme(travel~1,Rail,list(Rail=~1))

## 6.3.1 Tree growth: An example using `lme'

## NOTE: this example is extremely sensitive to changes in the 
## nlme package.

Loblolly$age <- Loblolly$age - mean(Loblolly$age)
lmc <- lmeControl(niterEM=500,msMaxIter=100)

m0 <- lme(height ~ age + I(age^2) + I(age^3),Loblolly,
          random=list(Seed=~age+I(age^2)+I(age^3)),
          correlation=corAR1(form=~age|Seed),control=lmc)
plot(m0)

m1 <- lme(height ~ age+I(age^2)+I(age^3)+I(age^4),Loblolly,
          list(Seed=~age+I(age^2)+I(age^3)),
          cor=corAR1(form=~age|Seed),control=lmc)
plot(m1)
m2 <- lme(height~age+I(age^2)+I(age^3)+I(age^4)+I(age^5),
          Loblolly,list(Seed=~age+I(age^2)+I(age^3)),
          cor=corAR1(form=~age|Seed),control=lmc)
plot(m2)

plot(m2,Seed~resid(.))
qqnorm(m2,~resid(.))
qqnorm(m2,~ranef(.))

m3 <- lme(height~age+I(age^2)+I(age^3)+I(age^4)+I(age^5),
          Loblolly,list(Seed=~age+I(age^2)+I(age^3)),control=lmc)
anova(m3,m2)

m4 <- lme(height~age+I(age^2)+I(age^3)+I(age^4)+I(age^5),
          Loblolly,list(Seed=~age+I(age^2)),
          correlation=corAR1(form=~age|Seed),control=lmc)
anova(m4,m2)


m5 <- lme(height~age+I(age^2)+I(age^3)+I(age^4)+I(age^5),
          Loblolly,list(Seed=pdDiag(~age+I(age^2)+I(age^3))),
          correlation=corAR1(form=~age|Seed),control=lmc)
anova(m2,m5)

plot(augPred(m2))

## 6.3.2 Several levels of nesting

lme(score~Machine,Machines,list(Worker=~1,Machine=~1))

## 6.5 GLMMs with R

library(MASS)

## NOTE: Have to run ch2.r analysis of these data first!

rf <- residuals(b4,type="d") # extract deviance residuals
## create an identifier for each sampling station
solr$station <- factor(with(solr,paste(-la,-lo,-t,sep="")))

## is there evidence of a station effect in the residuals?
solr$rf <-rf
rm <- lme(rf~1,solr,random=~1|station)
rm0 <- lm(rf~1,solr)
anova(rm,rm0)

b <- glmmPQL(eggs ~ offset(off)+lo+la+t+I(lo*la)+I(lo^2)+
            I(la^2)+I(t^2)+I(lo*t)+I(la*t)+I(lo^3)+I(la^3)+
            I(t^3)+I(lo*la*t)+I(lo^2*la)+I(lo*la^2)+I(lo^2*t)+
            I(la^2*t)+I(la*t^2)+I(lo*t^2) # end log spawn
            + a +I(a*t)+I(t^2*a),random=list(station=~1),
            family=quasi(link=log,variance="mu"),data=solr)
summary(b)

b1<-update(b,~.-I(lo*la*t))
summary(b1) 

b2<-update(b1,~.-I(lo*t))
summary(b2) 

b3<-update(b2,~.-I(lo^2*t))
summary(b3) 

b4<-update(b3,~.-I(la*t^2))
summary(b4) 

fv <- exp(fitted(b4)+solr$off) # note need to add offset
resid <- solr$egg-fv          # raw residuals
plot(fv^.5,solr$eggs^.5)
abline(0,1,lwd=2)
plot(fv^.5,resid/fv^.5)
plot(fv^.5,resid)
fl<-sort(fv^.5)
## add 1 s.d. and 2 s.d. reference lines
lines(fl,fl);lines(fl,-fl);lines(fl,2*fl,lty=2)
lines(fl,-2*fl,lty=2)

intervals(b4,which="var-cov")

## 6.7.1 A GAMM for sole eggs, note parametric `a' term removed,
## as s(t,k=5,by=a) no longer centred (since 1.4-0).

library(mgcv)

bam<-gamm(eggs~te(lo,la,t,bs=c("tp","tp"),k=c(25,5),d=c(2,1))
          +s(t,k=5,by=a)+offset(off),family=quasi(link=log,
          variance="mu"),data=solr,random=list(station=~1))

bam$gam

## 6.7.2 The temperature in Cairo

data(cairo)

ctamm<-gamm(temp~s(day.of.year,bs="cc",k=20)+s(time,bs="cr"),
         data=cairo,correlation=corAR1(form=~1|year))

summary(ctamm$gam)

intervals(ctamm$lme,which="var-cov")

ctamm$gam$sig2/ctamm$gam$sp

ctamm0 <- gamm(temp~s(day.of.year,bs="cc",k=20),data=cairo,
               correlation=corAR1(form=~1|year))
ctamm0$gam

anova(ctamm0$lme,ctamm$lme)
plot(ctamm$gam,scale=0)

