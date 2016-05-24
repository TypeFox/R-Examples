### R code from vignette source 'using-lsmeans.rnw'

###################################################
### code chunk number 1: using-lsmeans.rnw:66-68
###################################################
options(show.signif.stars=FALSE, prompt="R> ", continue="   ", 
    useFancyQuotes=FALSE, width=100, digits=6)


###################################################
### code chunk number 2: using-lsmeans.rnw:152-155
###################################################
library("lsmeans")
oranges.lm1 <- lm(sales1 ~ price1 + price2 + day + store, data = oranges)
anova(oranges.lm1)


###################################################
### code chunk number 3: using-lsmeans.rnw:158-159
###################################################
( oranges.rg1 <- ref.grid(oranges.lm1) )


###################################################
### code chunk number 4: using-lsmeans.rnw:162-163
###################################################
summary(oranges.rg1)


###################################################
### code chunk number 5: using-lsmeans.rnw:168-169
###################################################
lsmeans(oranges.rg1, "day")   ## or lsmeans(oranges.lm1, "day")


###################################################
### code chunk number 6: using-lsmeans.rnw:172-173
###################################################
with(oranges, tapply(sales1, day, mean))


###################################################
### code chunk number 7: using-lsmeans.rnw:184-186
###################################################
lsmeans(oranges.lm1, "day", at = list(price1 = 50, 
    price2 = c(40,60), day = c("2","3","4")) )


###################################################
### code chunk number 8: using-lsmeans.rnw:195-198
###################################################
org.lsm <- lsmeans(oranges.lm1, "day", by = "price2", 
    at = list(price1 = 50, price2 = c(40,60), day = c("2","3","4")) )
org.lsm


###################################################
### code chunk number 9: using-lsmeans.rnw:201-204 (eval = FALSE)
###################################################
## lsmeans(oranges.lm1, ~ day | price, at = ... )         # Ex 1
## lsmeans(oranges.lm1, c("day","price2"), at = ... )     # Ex 2
## lsmeans(oranges.lm1, ~ day * price, at = ... )         # Ex 3


###################################################
### code chunk number 10: using-lsmeans.rnw:214-215
###################################################
str(org.lsm)


###################################################
### code chunk number 11: using-lsmeans.rnw:221-223
###################################################
( org.sum <- summary(org.lsm, infer = c(TRUE,TRUE), 
                    level = .90, adjust = "bon", by = "day") )


###################################################
### code chunk number 12: using-lsmeans.rnw:228-229
###################################################
class(org.sum)


###################################################
### code chunk number 13: using-lsmeans.rnw:233-234
###################################################
transform(org.sum, lsrubles = lsmean * 34.2)


###################################################
### code chunk number 14: using-lsmeans.rnw:242-244
###################################################
org.lsm2 <- update(org.lsm, by.vars = NULL, level = .99)
org.lsm2


###################################################
### code chunk number 15: org-plot
###################################################
plot(org.lsm, by = "price2")


###################################################
### code chunk number 16: using-lsmeans.rnw:267-268
###################################################
contrast(org.lsm, method = "eff")


###################################################
### code chunk number 17: using-lsmeans.rnw:273-275
###################################################
days.lsm <- lsmeans(oranges.rg1, "day")
( days_contr.lsm <- contrast(days.lsm, "trt.vs.ctrl", ref = c(5,6)) )


###################################################
### code chunk number 18: using-lsmeans.rnw:280-281 (eval = FALSE)
###################################################
## confint(contrast(days.lsm, "trt.vs.ctrlk"))


###################################################
### code chunk number 19: using-lsmeans.rnw:289-290
###################################################
pairs(org.lsm)


###################################################
### code chunk number 20: using-lsmeans.rnw:293-294
###################################################
cld(days.lsm, alpha = .10)


###################################################
### code chunk number 21: days-cmp
###################################################
plot(days.lsm, comparisons = TRUE, alpha = .10)


###################################################
### code chunk number 22: using-lsmeans.rnw:321-322
###################################################
rbind(pairs(org.lsm)[1:3], pairs(org.lsm, by = "day")[1])


###################################################
### code chunk number 23: using-lsmeans.rnw:327-328
###################################################
rbind(pairs(lsmeans(org.lsm, "day")), pairs(lsmeans(org.lsm, "price2")))


###################################################
### code chunk number 24: using-lsmeans.rnw:334-337
###################################################
oranges.mlm <- lm(cbind(sales1,sales2) ~ price1 + price2 + day + store, 
                 data = oranges)
ref.grid(oranges.mlm)


###################################################
### code chunk number 25: using-lsmeans.rnw:340-342
###################################################
org.mlsm <- lsmeans(oranges.mlm, ~ day | variety, mult.name = "variety")
cld(org.mlsm, sort = FALSE)


###################################################
### code chunk number 26: using-lsmeans.rnw:350-351
###################################################
org.vardiff <- update(pairs(org.mlsm, by = "day"), by = NULL)


###################################################
### code chunk number 27: using-lsmeans.rnw:354-355
###################################################
cld(org.vardiff)


###################################################
### code chunk number 28: using-lsmeans.rnw:361-362
###################################################
contrast(org.mlsm, interaction = c("poly", "pairwise"))


###################################################
### code chunk number 29: using-lsmeans.rnw:370-372
###################################################
# Ensure we see the same results each time
set.seed(123454321)


###################################################
### code chunk number 30: using-lsmeans.rnw:374-377
###################################################
library("multcomp")
days.glht <- as.glht(days_contr.lsm)
summary(days.glht, test = adjusted("Westfall"))


###################################################
### code chunk number 31: using-lsmeans.rnw:380-382 (eval = FALSE)
###################################################
## days.glht1 <- glht(oranges.lm1, 
##                    lsm("day", contr = "trt.vs.ctrl", ref = c(5,6)))


###################################################
### code chunk number 32: using-lsmeans.rnw:386-388 (eval = FALSE)
###################################################
## summary(days_contr.lsm, adjust = "mvt")
## summary(days.glht)


###################################################
### code chunk number 33: using-lsmeans.rnw:394-395 (eval = FALSE)
###################################################
## summary(as.glht(pairs(org.lsm)))


###################################################
### code chunk number 34: using-lsmeans.rnw:398-399 (eval = FALSE)
###################################################
## summary(as.glht(pairs(org.lsm), by = NULL))


###################################################
### code chunk number 35: using-lsmeans.rnw:402-403 (eval = FALSE)
###################################################
## summary(as.glht(pairs(org.lsm, by = NULL)))


###################################################
### code chunk number 36: using-lsmeans.rnw:415-420
###################################################
data("Oats", package = "nlme")
library("lme4")
Oats.lmer <- lmer(log(yield) ~ Variety*factor(nitro) + (1|Block/Variety), 
                 data = Oats)
anova(Oats.lmer)


###################################################
### code chunk number 37: oatcontr (eval = FALSE)
###################################################
## contrast(lsmeans(Oats.lmer, "nitro"), "poly")


###################################################
### code chunk number 38: using-lsmeans.rnw:427-428
###################################################
cat("NOTE: Results may be misleading due to involvement in interactions")


###################################################
### code chunk number 39: using-lsmeans.rnw:430-431
###################################################
contrast(lsmeans(Oats.lmer, "nitro"), "poly")


###################################################
### code chunk number 40: using-lsmeans.rnw:435-437
###################################################
Oats.lmer2 <- lmer(log(yield) ~ Variety + poly(nitro,2) 
                                + (1|Block/Variety),  data = Oats)


###################################################
### code chunk number 41: using-lsmeans.rnw:441-442
###################################################
Oats.lsm2 <- lsmeans(Oats.lmer2, ~ nitro | Variety, cov.reduce = FALSE)


###################################################
### code chunk number 42: using-lsmeans.rnw:455-460
###################################################
library("xtable")
xtbl <- xtable(Oats.lsm2, caption = "Example using \\texttt{xtable}",
    label = "xtable:example")
print(xtbl, table.placement = "t")    
cat("See Table~\\ref{xtable:example}.\n")


###################################################
### code chunk number 43: oatslmer
###################################################
lsmip(Oats.lmer, Variety ~ nitro, ylab = "Observed log(yield)")


###################################################
### code chunk number 44: oatslmer2
###################################################
lsmip(Oats.lsm2, Variety ~ nitro, ylab = "Predicted log(yield)")


###################################################
### code chunk number 45: using-lsmeans.rnw:496-497
###################################################
str(Oats.lsm2)


###################################################
### code chunk number 46: using-lsmeans.rnw:500-501
###################################################
summary(Oats.lsm2, type = "response")


###################################################
### code chunk number 47: using-lsmeans.rnw:512-516
###################################################
Oats.log1 <- lmer(log(yield + 5) ~ Variety + factor(nitro) 
                  + (1|Block/Variety), data = Oats)
( Oats.rg1 <- update(ref.grid(Oats.log1), 
                    tran = make.tran("genlog", 5)) )


###################################################
### code chunk number 48: using-lsmeans.rnw:521-522
###################################################
round(predict(Oats.rg1, type = "response"), 1)


###################################################
### code chunk number 49: using-lsmeans.rnw:527-529
###################################################
my.tran <- make.tran("boxcox", c(.567, 10))
my.tran$linkfun(10:15)


###################################################
### code chunk number 50: using-lsmeans.rnw:537-541
###################################################
Oats.bc <- with(my.tran, lmer(linkfun(yield) ~ Variety + factor(nitro)
                              + (1|Block/Variety), data = Oats))
( rg.bc <- ref.grid(Oats.bc) )
round(predict(rg.bc, type = "response"), 1)


###################################################
### code chunk number 51: using-lsmeans.rnw:546-547
###################################################
rg.bc.regrid <- regrid(rg.bc)


###################################################
### code chunk number 52: using-lsmeans.rnw:550-551
###################################################
round(rg.bc.regrid@bhat, 1)


###################################################
### code chunk number 53: using-lsmeans.rnw:556-558
###################################################
summary(lsmeans(rg.bc, "Variety"), type = "response")
lsmeans(rg.bc.regrid, "Variety")


###################################################
### code chunk number 54: using-lsmeans.rnw:577-581
###################################################
rg.log <- regrid(rg.bc, "log")
lsm.log <- lsmeans(rg.log, "Variety")
summary(lsm.log, type = "response")
summary(pairs(lsm.log), type = "response")


###################################################
### code chunk number 55: using-lsmeans.rnw:592-594
###################################################
Oats.Vlsm = lsmeans(Oats.lmer2, "Variety")
test(Oats.Vlsm, null = log(100), type = "response")


###################################################
### code chunk number 56: using-lsmeans.rnw:606-607
###################################################
test(Oats.Vlsm, null = log(100), delta = 0.20, type = "r")


###################################################
### code chunk number 57: using-lsmeans.rnw:614-615
###################################################
test(contrast(Oats.Vlsm, "trt.vs.ctrlk"), side = ">")


###################################################
### code chunk number 58: using-lsmeans.rnw:619-620
###################################################
test(contrast(Oats.Vlsm, "trt.vs.ctrlk"), side = "nonsup", delta = .25)


###################################################
### code chunk number 59: chick-plot
###################################################
require("lattice")
xyplot(weight~Time | Diet, groups = ~ Chick, data = ChickWeight, 
    type = "o", layout=c(4, 1))


###################################################
### code chunk number 60: using-lsmeans.rnw:642-644
###################################################
Chick.lmer <- lmer(weight ~ Diet * Time + (0 + Time | Chick), 
    data = ChickWeight)


###################################################
### code chunk number 61: using-lsmeans.rnw:647-648
###################################################
( Chick.lst <- lstrends (Chick.lmer, ~ Diet, var = "Time") )


###################################################
### code chunk number 62: using-lsmeans.rnw:651-652
###################################################
cld (Chick.lst)


###################################################
### code chunk number 63: using-lsmeans.rnw:664-667
###################################################
lsm.options(ref.grid = list(level = .90),
            lsmeans = list(),
            contrast = list(infer = c(TRUE, TRUE)))


###################################################
### code chunk number 64: using-lsmeans.rnw:672-673
###################################################
get.lsm.option("estble.tol")


###################################################
### code chunk number 65: using-lsmeans.rnw:690-691
###################################################
lsmeans(Oats.lmer2, pairwise ~ Variety)


###################################################
### code chunk number 66: using-lsmeans.rnw:695-696
###################################################
lsm.options(ref.grid = NULL, contrast = NULL)


###################################################
### code chunk number 67: using-lsmeans.rnw:706-709
###################################################
nutr.lm <- lm(gain ~ (age + group + race)^2, data = nutrition)
library("car")
Anova(nutr.lm)


###################################################
### code chunk number 68: nutr-intplot
###################################################
lsmip(nutr.lm, race ~ age | group)
lsmeans(nutr.lm, ~ group*race)


###################################################
### code chunk number 69: using-lsmeans.rnw:726-728
###################################################
nutr.lsm <- lsmeans(nutr.lm, ~ group * race, weights = "proportional",
    at = list(age = c("2","3"), race = c("Black","White")))


###################################################
### code chunk number 70: using-lsmeans.rnw:731-734
###################################################
nutr.lsm    
summary(pairs(nutr.lsm, by = "race"), by = NULL)
summary(pairs(nutr.lsm, by = "group"), by = NULL)


###################################################
### code chunk number 71: using-lsmeans.rnw:747-751
###################################################
lsmeans(nutr.lm, "race", weights = "equal")
lsmeans(nutr.lm, "race", weights = "prop")
lsmeans(nutr.lm, "race", weights = "outer")
lsmeans(nutr.lm, "race", weights = "cells")


###################################################
### code chunk number 72: using-lsmeans.rnw:760-762
###################################################
temp = lsmeans(nutr.lm, c("group","race"), weights = "prop")
lsmeans(temp, "race", weights = "prop")


###################################################
### code chunk number 73: using-lsmeans.rnw:767-768
###################################################
with(nutrition, tapply(gain, race, mean))


###################################################
### code chunk number 74: using-lsmeans.rnw:776-780
###################################################
library("mediation")
levels(framing$educ) = c("NA","Ref","< HS", "HS", "> HS","Coll +")
framing.glm = glm(cong_mesg ~ age + income + educ + emo + gender * factor(treat),
                  family = binomial, data = framing)


###################################################
### code chunk number 75: framinga
###################################################
lsmip(framing.glm, treat ~ educ | gender, type = "response")


###################################################
### code chunk number 76: framingb
###################################################
lsmip(framing.glm, treat ~ educ | gender, type = "response",
      cov.reduce = emo ~ treat*gender + age + educ + income)


###################################################
### code chunk number 77: using-lsmeans.rnw:808-810
###################################################
ref.grid(framing.glm, 
    cov.reduce = emo ~ treat*gender + age + educ + income)@grid


###################################################
### code chunk number 78: using-lsmeans.rnw:841-843 (eval = FALSE)
###################################################
## rg <- ref.grid(my.model, at = list(x1 = c(5,10,15)),
##                cov.reduce = list(x2 ~ x1,  x3 ~ x1 + x2))


###################################################
### code chunk number 79: housing-plot
###################################################
library("ordinal")
data(housing, package = "MASS")
housing.clm <- clm(Sat ~ (Infl + Type + Cont)^2,
                   data = housing, weights = Freq, link = "probit")
lsmip(housing.clm, Cont ~ Infl | Type, layout = c(4,1))


###################################################
### code chunk number 80: using-lsmeans.rnw:895-896
###################################################
test(pairs(lsmeans(housing.clm, ~ Infl | Type)), joint = TRUE)


###################################################
### code chunk number 81: using-lsmeans.rnw:899-900
###################################################
test(pairs(lsmeans(housing.clm, ~ Cont | Type)), joint = TRUE)


###################################################
### code chunk number 82: using-lsmeans.rnw:905-906
###################################################
ref.grid(housing.clm, mode = "cum.prob")


###################################################
### code chunk number 83: using-lsmeans.rnw:909-911
###################################################
lsmeans(housing.clm, ~ Infl, at = list(cut = "Medium|High"), 
        mode = "cum.prob")


###################################################
### code chunk number 84: using-lsmeans.rnw:914-916
###################################################
summary(lsmeans(housing.clm, ~ Infl, at = list(cut = "Medium|High"), 
                mode = "linear.predictor"), type = "response")


###################################################
### code chunk number 85: using-lsmeans.rnw:924-932
###################################################
require("nlme")
options(contrasts = c("contr.treatment", "contr.poly"))
Chick.nlme = nlme(weight ~ SSlogis(Time, asym, xmid, scal), 
    data = ChickWeight,
    fixed = list(asym + xmid ~ Diet, scal ~ 1),
    random = asym ~ 1 | Chick, 
    start = c(200, 100, 200, 100,   10, 0, 0, 0,   7))
Chick.nlme


###################################################
### code chunk number 86: using-lsmeans.rnw:935-937
###################################################
cld(lsmeans(Chick.nlme, ~ Diet, param = "asym"))    
cld(lsmeans(Chick.nlme, ~ Diet, param = "xmid"))    


###################################################
### code chunk number 87: using-lsmeans.rnw:949-954
###################################################
library("MCMCpack")
counts <- c(18, 17, 15,   20, 10, 20,   25, 13, 12)
outcome <- gl(3, 1, 9)
treatment <- gl(3, 3)
posterior <- MCMCpoisson(counts ~ outcome + treatment, mcmc = 1000)


###################################################
### code chunk number 88: using-lsmeans.rnw:957-958
###################################################
( post.lsm <- lsmeans(posterior, "treatment") )


###################################################
### code chunk number 89: using-lsmeans.rnw:961-963
###################################################
library("coda")
summary(as.mcmc(post.lsm))


