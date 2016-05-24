names(eba1977)
attach(eba1977)
fit <- glm(cases~city+age+offset(log(pop)), family=poisson)
summary(fit)
min(fitted(fit))
pchisq(deviance(fit), df.residual(fit), lower=F)
pchisq(23.45, 15, lower=F)
drop1(fit, test="Chisq")
fit2 <- glm(cases~(city=="Fredericia")+age+offset(log(pop)), 
                 family=poisson)
anova(fit, fit2, test="Chisq")
drop1(fit2,  test="Chisq")
summary(fit2)
cf <- coefficients(summary(fit2))
est <- cf[,1]
s.e. <- cf[,2]
rr <- exp(cbind(est, est - s.e.*qnorm(.975), est 
                     + s.e.*qnorm(.975) ))
colnames(rr) <- c("RateRatio", "CI.lo","CI.hi")
rr
exp(cbind(coef(fit2), confint(fit2)))
head(nickel.expand)
subset(nickel.expand, id==325)
nickel.expand <- within(nickel.expand, 
     lung.cancer <- as.numeric(icd %in% c(162,163)))
attach(nickel.expand)
pyr <- tapply(ageout-agein,list(ygr,agr), sum)
print(round(pyr), na.print="-")
count <- tapply(lung.cancer, list(ygr, agr), sum)
print(count, na.print="-")
print(round(count/pyr*1000, 1), na.print="-")
expect.count <- tapply(lung/1e6*(ageout-agein),
                       list(ygr,agr), sum)
print(round(expect.count, 1), na.print="-")
expect.tot <- sum(lung/1e6*(ageout-agein))
expect.tot
count.tot <- sum(lung.cancer)
count.tot
count.tot/expect.tot
fit <- glm(lung.cancer ~ 1, poisson,
           offset = log((ageout-agein)*lung/1e6))
summary(fit)
exp(coef(fit))
tapply(lung.cancer, agr, sum)
tapply(lung.cancer, ygr, sum)
detach()
nickel.expand <- within(nickel.expand,{
    A <- factor(agr)
    Y <- factor(ygr)
    lv <- levels(A)
    lv[1:6] <- "< 50"
    lv[11:13] <- "70+"
    levels(A) <- lv
    lv <- levels(Y)
    lv[7:10] <- "1961ff"
    levels(Y) <- lv
    rm(lv)
})
attach(nickel.expand)
fit <- glm(lung.cancer ~ A + Y, poisson,
           offset=log((ageout-agein)*lung/1e6))
drop1(fit, test="Chisq")
fit <- glm(lung.cancer ~ Y - 1, poisson,
           offset=log((ageout-agein)*lung/1e6))
summary(fit)
round(exp(coef(fit)), 1)
expect.count <-  tapply(lung/1e6*(ageout-agein), Y, sum)
count <- tapply(lung.cancer, Y, sum)
cbind(count=count, expect=round(expect.count,1),
      SMR= round(count/expect.count, 1))
detach()
nickel.expand <- within(nickel.expand,{
    TFE <- cut(agein-age1st, c(0,20,30,40,50,100), right=F)
    AFE <- cut(age1st, c(0, 20, 27.5, 35, 100), right=F)
    YFE <- cut(dob + age1st, c(0, 1910, 1915, 1920, 1925),right=F)
    EXP <- cut(exposure, c(0, 0.5, 4.5, 8.5, 12.5, 25), right=F)
})
attach(nickel.expand)
fit <- glm(lung.cancer ~ TFE + AFE + YFE + EXP, poisson,
           offset=log((ageout-agein)*lung/1e6))
drop1(fit, test="Chisq")
summary(fit)
rm(list=ls())
while(search()[2] != "package:ISwR") detach()
