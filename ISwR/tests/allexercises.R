library(ISwR)
ps.options(height=3.5, width=4.4, pointsize=8, horiz=F)
par(mar=c(4,4,3,2)+.1)
set.seed(310367)
x <- y <- c(7, 9, NA, NA, 13)
all(is.na(x) == is.na(y)) & all((x == y)[!is.na(x)])
x <- factor(c("Huey", "Dewey", "Louie", "Huey"))
y <- c("blue", "red", "green")
x
y[x]
juul.girl <- juul[juul$age >=7 & juul$age < 14 & juul$sex == 2,]
summary(juul.girl)
sapply(1:10, function(i) mean(rexp(20)))
rm(list=ls())
while(search()[2] != "package:ISwR") detach()
x <- 1:10
z <- append(x, 1.23, after=7)
z
z <- c(x[1:7],1.23,x[8:10])
z
v <- 1.23; k <- 7
i <- seq(along=x)
z <- c(x[i <= k], v, x[i > k])
z
write.table(thuesen, file="foo.txt")
# edit the file
read.table("foo.txt", na.strings=".")
write.table(thuesen, file="foo.txt", na=".")
read.table("foo.txt", na.strings=".")
rm(list=ls())
while(search()[2] != "package:ISwR") detach()
1 - pnorm(3)
1 - pnorm(42, mean=35, sd=6)
dbinom(10, size=10, prob=0.8)
punif(0.9) # this one is obvious...
1 - pchisq(6.5, df=2)
pnorm(-2) * 2
qnorm(1-.01/2)
qnorm(1-.005/2)
qnorm(1-.001/2)
qnorm(.25)
qnorm(.75)
rbinom(10, 1, .5)
ifelse(rbinom(10, 1, .5) == 1, "H", "T")
c("H", "T")[1 + rbinom(10, 1, .5)]
rm(list=ls())
while(search()[2] != "package:ISwR") detach()
x <- 1:5 ; y <- rexp(5,1) ; opar <- par(mfrow=c(2,2))
plot(x, y, pch=15) # filled square
plot(x, y, type="b", lty="dotted")
plot(x, y, type="b", lwd=3)
plot(x, y, type="o", col="blue")
par(opar)
plot(rnorm(10),type="o", pch=21, bg="white")
x1 <- rnorm(20)
x2 <- rnorm(10)+1
q1 <- qqnorm(x1, plot.it=F)
q2 <- qqnorm(x2, plot.it=F)
xr <- range(q1$x, q2$x)
yr <- range(q1$y, q2$y)
qqnorm(x1, xlim=xr, ylim=yr)
points(q2, col="red")
hist(react)
library(MASS)
truehist(react,h=1,x0=.5)
z <- runif(5)
curve(quantile(z,x), from=0, to=1)
rm(list=ls())
while(search()[2] != "package:ISwR") detach()
qqnorm(react)
t.test(react)
wilcox.test(react)
wilcox.test(vital.capacity~group, data=vitcap)
attach(intake) ; opar <- par(mfrow=c(2,2))
plot(post ~ pre) ; abline(0,1)
plot((post+pre)/2, post - pre,
    ylim=range(0,post-pre)); abline(h=0)
hist(post-pre)
qqnorm(post-pre)
detach(intake)
par(opar)
shapiro.test(react)
shapiro.test(react[-c(1,334)])
qqnorm(react[-c(1,334)])
attach(ashina)
t.test(vas.active, vas.plac, paired=TRUE)
t.test((vas.active-vas.plac)[grp==1],
       (vas.plac-vas.active)[grp==2])
t.test(rnorm(25))$p.value       #repeat 10x
t.test(rt(25,df=2))$p.value     #repeat 10x
t.test(rexp(25), mu=1)$p.value  #repeat 10x
x <- replicate(5000, t.test(rexp(25), mu=1)$p.value)
qqplot(sort(x),ppoints(5000),type="l",log="xy")
rm(list=ls())
while(search()[2] != "package:ISwR") detach()
fit <- lm(metabolic.rate ~ body.weight, data=rmr)
summary(fit)
811.2267 + 7.0595 * 70 # , or:
predict(fit, newdata=data.frame(body.weight=70))
qt(.975,42)
7.0595 + c(-1,1) * 2.018 * 0.9776 # , or:
confint(fit)
summary(lm(log(ab)~age, data=malaria))
plot(log(ab)~age, data=malaria)
rho <- .90 ; n <- 100
x <- rnorm(n)
y <- rnorm(n, rho * x, sqrt(1 - rho^2))
plot(x, y)
cor.test(x, y)
cor.test(x, y, method="spearman")
cor.test(x, y, method="kendall")
rm(list=ls())
while(search()[2] != "package:ISwR") detach()
walk <- unlist(zelazo) # or c(..,recursive=TRUE)
group <- factor(rep(1:4,c(6,6,6,5)), labels=names(zelazo))
summary(lm(walk ~ group))
t.test(zelazo$active,zelazo$ctr.8w) # first vs. last
t.test(zelazo$active,unlist(zelazo[-1])) # first vs. rest
fit <- lm(volume~method+subject, data=lung)
anova(fit)
summary(fit)
kruskal.test(walk ~ group)
wilcox.test(zelazo$active,zelazo$ctr.8w) # first vs. last
wilcox.test(zelazo$active,unlist(zelazo[-1])) # first vs. rest
friedman.test(volume ~ method | subject, data=lung)
wilcox.test(lung$volume[lung$method=="A"],
            lung$volume[lung$method=="C"], paired=TRUE) # etc.
attach(juul)
tapply(sqrt(igf1),tanner, sd, na.rm=TRUE)
plot(sqrt(igf1)~jitter(tanner))
oneway.test(sqrt(igf1)~tanner)
rm(list=ls())
while(search()[2] != "package:ISwR") detach()
binom.test(0, 10, p=.20, alt="less")
binom.test(0, 13, p=.20, alt="less")
binom.test(0, 14, p=.20, alt="less")
prop.test(c(210,122),c(747,661))
M <- matrix(c(23,7,18,13),2,2)
chisq.test(M)
fisher.test(M)
prop.test(M)
tbl <- c(42, 157, 47, 62, 4, 15, 4, 1, 8, 28, 9, 7)
dim(tbl) <- c(2,2,3)
dimnames(tbl) <- list(c("A","B"),
                      c("not pierced","pierced"),
                      c("ok","broken","cracked"))
ftable(tbl)
fisher.test(tbl["B",,]) # slice analysis
fisher.test(tbl["A",,])
fisher.test(margin.table(tbl,2:3)) # marginal
p <- seq(0,1,0.001)
pval <- sapply(p,function(p)binom.test(3,15,p=p)$p.value)
plot(p,pval,type="l")
rm(list=ls())
while(search()[2] != "package:ISwR") detach()
power.t.test(power=.8,delta=.30,sd=.20)
power.t.test(power=.8,delta=.30,sd=.20,alt="one.sided")
(qnorm(.975)+qnorm(.8))^2*2*(.2/.3)^2 # approx. formula
power.t.test(n=8, delta=.30, sd=.20)   # power with eq.size
d2 <- .30 * sqrt(2/8) / sqrt(1/6+1/10) # corr.f.uneq. size
power.t.test(n=8, delta=d2, sd=.20)
power.prop.test(power=.9, p1=.6, p2=.75)
power.prop.test(power=.8, p1=.6, p2=.75)
curve(dt(x-3, 25), from=0, to=5)
curve(dt(x, 25, 3), add=TRUE)
rm(list=ls())
while(search()[2] != "package:ISwR") detach()
attach(thuesen)
f <- cut(blood.glucose, c(4, 7, 9, 12, 20))
levels(f) <- c("low", "intermediate", "high", "very high")
bcmort2 <- within(bcmort,{
  period <- area <- cohort
  levels(period) <- rep(c("1991-2001","1981-1991"), each=2)
  levels(area) <- rep(c("Cph+Frb","Nat"),2)
})
summary(bcmort2)
ashina.long <- reshape(ashina, direction="long",
                    varying=1:2, timevar="treat")
ashina.long <- within(ashina.long, {
     m <- matrix(c(2,1,1,2),2)
     id <- factor(id)
     treat <- factor(treat)
     grp <- factor(grp)
     period <- factor(m[cbind(grp,treat)])
     rm(m)
})
within(ashina.long,
  period2 <- ifelse(treat != "active",
             as.numeric(grp), 3 - as.numeric(grp))
)
stroke.trim <- function(t1, t2)
   subset(transform(stroke,
                    entry=t1, exit=pmin(t2, obsmonths),
                    dead=dead & obsmonths <= t2),
          entry < exit)
stroke2 <- do.call(rbind, mapply(stroke.trim,
       c(0,0.5,2,12), c(0.5,2,12,Inf), SIMPLIFY=F))
table(stroke$dead)
table(stroke2$dead)
rm(list=ls())
while(search()[2] != "package:ISwR") detach()
summary(lm(log(bwt) ~ log(bpd) + log(ad), data=secher))
summary(lm(log(bwt) ~ log(ad), data=secher))
pairs(tlc)
summary(lm(log(tlc) ~ ., data=tlc))
opar <- par(mfrow=c(2,2))
plot(lm(log(tlc) ~ ., data=tlc), which=1:4)

drop1(lm(log(tlc) ~ ., data=tlc))
drop1(lm(log(tlc) ~ . - age, data=tlc))

par(mfrow=c(1,1))
plot(log(tlc) ~ height, data=tlc)
par(mfrow=c(2,2))
plot(lm(tlc ~ ., data=tlc), which=1:4) # slightly worse
par(opar)
summary(lm(sqrt(igf1) ~ age, data=juul2, subset=(age >= 25)))
anova(lm(sqrt(igf1) ~ age + weight + height,
           data=juul2, subset=(age >= 25)))
summary(lm(dl.milk ~ . - no, data=kfm))
summary(lm(dl.milk ~ . - no - mat.weight, data=kfm))
summary(lm(dl.milk ~ . - no - mat.weight - sex, data=kfm))
summary(lm(dl.milk ~ weight + mat.height, data=kfm))
rm(list=ls())
while(search()[2] != "package:ISwR") detach()
ashina.long <- reshape(ashina, direction="long",
                    varying=1:2, timevar="treat")
ashina.long <- within(ashina.long, {
     m <- matrix(c(2,1,1,2),2)
     id <- factor(id)
     treat <- factor(treat)
     grp <- factor(grp)
     period <- factor(m[cbind(grp,treat)])
     rm(m)
})
fit.ashina <- lm(vas ~ id + period + treat, data=ashina.long)
drop1(fit.ashina, test="F")
anova(fit.ashina)

attach(ashina)
dd <- vas.active - vas.plac
t.test(dd[grp==1], -dd[grp==2], var.eq=T)
t.test(dd[grp==1], dd[grp==2], var.eq=T)
attach(tb.dilute)
anova(lm(reaction ~ animal + logdose))
ld <- c(0.5, 0, -0.5)[logdose]
anova(lm(reaction ~ animal + ld))
summary(lm(reaction ~ animal + ld))
4.7917 + 0.6039 * qt(c(.025,.975), 11)
# or:
confint(lm(reaction ~ animal + ld))["ld",]

slopes <- reaction[logdose==0.5] - reaction[logdose==-0.5]
t.test(slopes)

anova(lm(reaction ~ animal*ld))
a <- gl(2, 2, 8)
b <- gl(2, 4, 8)
x <- 1:8
y <- c(1:4,8:5)
z <- rnorm(8)
model.matrix(~ a:b)         ; lm(z ~ a:b)
model.matrix(~ a * b)       ; lm(z ~ a * b)
model.matrix(~ a:x)         ; lm(z ~ a:x)
model.matrix(~ a * x)       ; lm(z ~ a * x)
model.matrix(~ b * (x + y)) ; lm(z ~ b * (x + y))
attach(secretin)
model1 <- lm(gluc ~ person * time)
model2 <- lm(gluc ~ person + time)
model3 <- lm(gluc ~ person * time20plus + time)
model4 <- lm(gluc ~ person * time20plus + time.comb)
tt <- c(20,30,60,90,0)[time]
plot(fitted(model4)~tt,pch=as.character(person))
bp.obese <- transform(bp.obese,sex=factor(sex, labels=c("M","F")))
plot(log(bp) ~ log(obese), pch=c(20,21)[sex], data=bp.obese)
summary(lm(log(bp) ~ sex, data=bp.obese))
summary(lm(log(bp) ~ sex + log(obese), data=bp.obese))
summary(lm(log(bp) ~ sex*log(obese), data=bp.obese))
vitcap2 <- transform(vitcap2,group=factor(group,
                                          labels=c("exp>10",
                                          "exp<10", "unexp")))
attach(vitcap2)
plot(vital.capacity~age, pch=(20:22)[group])
vit.fit <- lm(vital.capacity ~ age*group)
summary(vit.fit)
drop1(vit.fit, test="F")
for (i in 1:3) abline(lm(vital.capacity ~ age,
                         subset=as.numeric(group)==i), lty=i)
legend(20, 3.5 ,legend=levels(group), pch=20:22, lty=1:3)
juul.prepub <- subset(juul, tanner==1)

summary(lm(sqrt(igf1)~age, data=juul.prepub, subset= sex==1))
summary(lm(sqrt(igf1)~age, data=juul.prepub, subset= sex==2))

summary(lm(sqrt(igf1)~age*factor(sex), data=juul.prepub))
summary(lm(sqrt(igf1)~age+factor(sex), data=juul.prepub))
summary(fit.aicopt <- step(lm(dl.milk ~ . - no, data=kfm)))
opar <- par(mfrow=c(2,2))
plot(fit.aicopt, which=1:4)
kfm[32,]
summary(kfm)
summary(update(fit.aicopt, ~ . - sex))
plot(update(fit.aicopt, ~ . - sex - ml.suppl), which=1:4)
par(opar)
juulyoung <- subset(juul, age < 25)
juulyoung <- transform(juulyoung,
                 sex=factor(sex), tanner=factor(tanner))
fit.untf <- lm(igf1 ~ age * sex * tanner, data=juulyoung,
               na.action=na.exclude)
plot(fitted(fit.untf) ~ age, data=juulyoung,
     col=c("red","green")[sex])
fit.log <- update(fit.untf, log(igf1) ~ .)
fit.sqrt <- update(fit.untf, sqrt(igf1) ~ .)
opar <- par(mfrow=c(2,2))
plot(fit.untf, which=1:4)
plot(fit.log, which=1:4)
plot(fit.sqrt, which=1:4)
par(opar)
rm(list=ls())
while(search()[2] != "package:ISwR") detach()
summary(glm(mal~age+log(ab), binomial, data=malaria))
attach(graft.vs.host)
type <- factor(type,labels=c("AML", "ALL", "CML"))
m1 <- glm(gvhd~rcpage+donage+type+preg+log(index), binomial)
m1a <- glm(gvhd~rcpage+donage+type+preg+index, binomial)
summary(m1)
summary(m1a)
drop1(m1, test="Chisq")
drop1(update(m1, ~ . - rcpage), test="Chisq")
drop1(update(m1, ~ . - rcpage - type), test="Chisq")
drop1(update(m1, ~ . - rcpage - type - preg), test="Chisq")
summary(m2 <- glm(gvhd~donage + log(index), binomial))
confint(m2)
## normal approximation:
est <-  coefficients(summary(m2))[,1]
se <-  coefficients(summary(m2))[,2]
est + cbind(qnorm(.025)*se, qnorm(.975)*se)
confint.default(m2)
counts <- c(13,40,157,40,21,61)
total <- c(108,264,375,310,181,162)
age <- gl(3,1,6)
type <- gl(2,3,6)
anova(glm(counts/total~age+type,weights=total, binomial),
      test="Chisq")
juul.girl <- transform(subset(juul,age>8 & age<20 &
                               complete.cases(menarche)),
                       menarche=factor(menarche))
logit.menarche <- glm(menarche~age+I(age^2)+I(age^3),
                      binomial, data=juul.girl)
probit.menarche <- glm(menarche~age+I(age^2)+I(age^3),
                       binomial(probit), data=juul.girl)
summary(logit.menarche)
summary(probit.menarche)
Age=seq(8,20,.1)
newages <- data.frame(age=Age)
p.logit <- predict(logit.menarche,newdata=newages,type="resp")
p.probit <- predict(probit.menarche,newdata=newages,type="resp")
matplot(Age,cbind(p.probit,p.logit),type="l")
rm(list=ls())
while(search()[2] != "package:ISwR") detach()
library(survival)
attach(graft.vs.host)
plot(survfit(Surv(time,dead)~gvhd))
survdiff(Surv(time,dead)~gvhd)
summary(coxph(Surv(time,dead) ~ gvhd)) # for comparison
summary(coxph(Surv(time,dead) ~
              gvhd + log(index) + donage + rcpage + preg))
attach(melanom)
cox1 <- coxph(Surv(days, status==1) ~
              log(thick) + sex + strata(ulc))
new <- data.frame(sex=2, thick=c(0.1, 0.2, 0.5))
svfit <-  survfit(cox1,newdata=new)
plot(svfit[2], ylim=c(.985, 1))
summary(coxph(Surv(obsmonths, dead)~age+sex, data=stroke))
summary(coxph(Surv(obsmonths, dead)~sex, data=stroke))
with(stroke, tapply(age,sex,mean))
stroke.trim <- function(t1, t2)
   subset(transform(stroke,
                    entry=t1, exit=pmin(t2, obsmonths),
                    dead=dead & obsmonths <= t2),
          entry < exit)
stroke2 <- do.call(rbind, mapply(stroke.trim,
       c(0,0.5,2,12), c(0.5,2,12,Inf), SIMPLIFY=F))
summary(coxph(Surv(entry, exit, dead)~age+sex, data=stroke2))
rm(list=ls())
while(search()[2] != "package:ISwR") detach()
bcmort2 <- within(bcmort,{
  period <- area <- cohort
  levels(period) <- rep(c("1991-2001","1981-1991"), each=2)
  levels(area) <- rep(c("Cph+Frb","Nat"),2)
})
bcfit <- glm(bc.deaths ~ (age + period + area)^2, poisson,
             offset=log(p.yr), data=bcmort2)
summary(bcfit)
drop1(bcfit, test="Chisq")
confint(bcfit, parm="period1981-1991:areaNat")
stroke.trim <- function(t1, t2)
   subset(transform(stroke,
                    entry=t1, exit=pmin(t2, obsmonths),
                    dead=dead & obsmonths <= t2),
          entry < exit)
stroke2 <- do.call(rbind, mapply(stroke.trim,
       c(0,0.5,2,12), c(0.5,2,12,Inf), SIMPLIFY=F))
summary(glm(dead~sex+age+factor(entry), poisson,
       offset=log(exit-entry), data=stroke2))
rm(list=ls())
while(search()[2] != "package:ISwR") detach()
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
# p <- profile(fit, alphamax=.2)
# alphamax semantics got changed in R 2.8
p <- profile(fit)
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
