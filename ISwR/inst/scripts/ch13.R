no.yes <- c("No","Yes")
smoking <- gl(2,1,8,no.yes)
obesity <- gl(2,2,8,no.yes)
snoring <- gl(2,4,8,no.yes)
n.tot <- c(60,17,8,2,187,85,51,23)
n.hyp <- c(5,2,1,0,35,13,15,8)
data.frame(smoking,obesity,snoring,n.tot,n.hyp)
expand.grid(smoking=no.yes, obesity=no.yes, snoring=no.yes)
hyp.tbl <- cbind(n.hyp,n.tot-n.hyp)
hyp.tbl 
glm(hyp.tbl~smoking+obesity+snoring,family=binomial("logit"))
glm(hyp.tbl~smoking+obesity+snoring,binomial)
prop.hyp <- n.hyp/n.tot
glm.hyp <- glm(prop.hyp~smoking+obesity+snoring,
               binomial,weights=n.tot)
glm(hyp.tbl~smoking+obesity+snoring, binomial("logit"))
glm.hyp <- glm(hyp.tbl~smoking+obesity+snoring,binomial)
summary(glm.hyp)
summary(glm(formula = hyp.tbl ~ smoking + obesity + snoring, family =
binomial))
glm.hyp <- glm(hyp.tbl~obesity+snoring,binomial) 
summary(glm.hyp)
glm.hyp <- glm(hyp.tbl~smoking+obesity+snoring,binomial)
anova(glm.hyp, test="Chisq")
glm.hyp <- glm(hyp.tbl~snoring+obesity+smoking,binomial)
anova(glm.hyp, test="Chisq")
glm.hyp <- glm(hyp.tbl~obesity+snoring,binomial)
anova(glm.hyp, test="Chisq")
drop1(glm.hyp, test="Chisq")
caesar.shoe
shoe.score <- 1:6                             
shoe.score
summary(glm(t(caesar.shoe)~shoe.score,binomial))
anova(glm(t(caesar.shoe)~shoe.score,binomial))
caesar.shoe.yes <- caesar.shoe["Yes",]
caesar.shoe.no <- caesar.shoe["No",]
caesar.shoe.total <- caesar.shoe.yes+caesar.shoe.no
prop.trend.test(caesar.shoe.yes,caesar.shoe.total)
prop.test(caesar.shoe.yes,caesar.shoe.total)
confint(glm.hyp)
confint.default(glm.hyp)
library(MASS)
plot(profile(glm.hyp))
if (.make.epsf) dev.copy2eps(file="profile-hyp.ps")
exp(cbind(OR=coef(glm.hyp), confint(glm.hyp)))
juul$menarche <- factor(juul$menarche, labels=c("No","Yes"))
juul$tanner <- factor(juul$tanner)
juul.girl <- subset(juul,age>8 & age<20 & 
                    complete.cases(menarche))
attach(juul.girl)
summary(glm(menarche~age,binomial))
summary(glm(menarche~age+tanner,binomial))
drop1(glm(menarche~age+tanner,binomial),test="Chisq")
predict(glm.hyp)
predict(glm.hyp, type="response")
glm.menarche <- glm(menarche~age, binomial)
Age <- seq(8,20,.1)
newages <- data.frame(age=Age)
predicted.probability <- predict(glm.menarche,
                                 newages,type="resp")
plot(predicted.probability ~ Age, type="l")
if (.make.epsf) dev.copy2eps(file="menarche-fit.ps")
fitted(glm.hyp)
prop.hyp   
fitted(glm.hyp)*n.tot
data.frame(fit=fitted(glm.hyp)*n.tot,n.hyp,n.tot)
age.group <- cut(age,c(8,10,12,13,14,15,16,18,20))
tb <- table(age.group,menarche)
tb
rel.freq <- prop.table(tb,1)[,2]
rel.freq
points(rel.freq ~ c(9,11,12.5,13.5,14.5,15.5,17,19),pch=5)
if (.make.epsf) dev.copy2eps(file="menarche-fit+obs.ps")
age.gr <- cut(age,c(8,12,13,14,20))
summary(glm(menarche~age+age.gr,binomial))
anova(glm(menarche~age+age.gr,binomial))
1-pchisq(8.058,3)
anova(glm(menarche~age+I(age^2)+I(age^3)+age.gr,binomial))
glm.menarche <- glm(menarche~age+I(age^2)+I(age^3), binomial)
predicted.probability <- 
    predict(glm.menarche, newages, type="resp")
plot(predicted.probability ~ Age, type="l")
points(rel.freq~c(9,11,12.5,13.5,14.5,15.5,17,19), pch=5)
if (.make.epsf) dev.copy2eps(file="menarche-cubic.ps")
rm(list=ls())
while(search()[2] != "package:ISwR") detach()
