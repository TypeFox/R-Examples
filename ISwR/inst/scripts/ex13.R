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
