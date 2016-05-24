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
