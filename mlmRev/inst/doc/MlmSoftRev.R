### R code from vignette source 'MlmSoftRev.Rnw'

###################################################
### code chunk number 1: preliminaries
###################################################
options(width=80, show.signif.stars = FALSE,
        lattice.theme = function() canonical.theme("pdf", color = FALSE))
library(mlmRev)
library(lme4)
library(lattice)
set.seed(1234321)


###################################################
### code chunk number 2: Examprep
###################################################
lmer(normexam ~ standLRT + sex + schgend + (1|school), Exam)


###################################################
### code chunk number 3: ExamData
###################################################
str(Exam)
summary(Exam)


###################################################
### code chunk number 4: ExamFit
###################################################
(Em1 <- lmer(normexam ~ standLRT + sex + schgend + (1|school), Exam))


###################################################
### code chunk number 5: Examtime
###################################################
system.time(lmer(normexam ~ standLRT + sex + schgend + (1|school), Exam))


###################################################
### code chunk number 6: ExamRelevel
###################################################
Exam$sex     <- relevel(Exam$sex, "M")
Exam$schgend <- relevel(Exam$schgend, "girls")
(Em2 <- lmer(normexam ~ standLRT + sex + schgend + (1|school), Exam))


###################################################
### code chunk number 7: ExamIds
###################################################
Exam <- within(Exam, ids <- factor(school:student))
str(Exam)


###################################################
### code chunk number 8: dupExamIds
###################################################
as.character(Exam$ids[which(duplicated(Exam$ids))])


###################################################
### code chunk number 9: OnlyBoy
###################################################
subset(Exam, ids == '43:86')
xtabs(~ sex + school, Exam, subset = school %in% c(43, 50, 52), drop = TRUE)


###################################################
### code chunk number 10: ExamXtabs
###################################################
xtabs(~ sex + school, Exam, subset = type == "Mxd", drop = TRUE)


###################################################
### code chunk number 11: ExamMosaicshow (eval = FALSE)
###################################################
## ExamMxd <- within(subset(Exam, type == "Mxd"), school <- factor(school))
## mosaicplot(~ school + sex, ExamMxd)


###################################################
### code chunk number 12: ExamMosaic
###################################################
ExamMxd <- within(subset(Exam, type == "Mxd"), school <- factor(school))
mosaicplot(~ school + sex, ExamMxd)


###################################################
### code chunk number 13: Examplot1
###################################################
print(xyplot(normexam ~ standLRT | sex * type, Exam,
             type = c("g", "p", "smooth"), layout = c(2,2),
             xlab = "Standardized London Reading Test score",
             ylab = "Normalized exam score", aspect = 1.2))


###################################################
### code chunk number 14: Examplot1show (eval = FALSE)
###################################################
## xyplot(normexam ~ standLRT | sex * type, Exam, type = c("g", "p", "smooth"))


###################################################
### code chunk number 15: Examplot2
###################################################
print(xyplot(normexam ~ standLRT, Exam, groups = sex:type,
             type = c("g", "smooth"), xlim = c(-3,3), ylim = c(-2,2),
             xlab = "Standardized London Reading Test score",
             ylab = "Normalized exam score",
             auto.key = list(space = 'right', points = FALSE, lines = TRUE),
             aspect = 1))


###################################################
### code chunk number 16: Examplot2show (eval = FALSE)
###################################################
## xyplot(normexam ~ standLRT, Exam, groups = sex:type, type = c("g", "smooth"))


###################################################
### code chunk number 17: Examplot3
###################################################
print(xyplot(normexam ~ standLRT | school, Exam,
             type = c("g", "p", "r"),
             xlab = "Standardized London Reading Test score",
             ylab = "Normalized exam score",
             subset = sex == "F" & type == "Sngl"))


###################################################
### code chunk number 18: Examplot3show
###################################################
xyplot(normexam ~ standLRT | school, Exam,
             type = c("g", "p", "r"),
             subset = sex == "F" & type == "Sngl")


###################################################
### code chunk number 19: Examplot4show (eval = FALSE)
###################################################
## xyplot(normexam ~ standLRT | school, Exam, type = c("g", "p", "r"),
##              subset = sex == "F" & type == "Sngl" & school != 48,
##              index.cond = function(x, y) coef(lm(y ~ x))[1])


###################################################
### code chunk number 20: Examplot4
###################################################
print(xyplot(normexam ~ standLRT | school, Exam,
             type = c("g", "p", "r"),
             xlab = "Standardized London Reading Test score",
             ylab = "Normalized exam score",
             subset = sex == "F" & type == "Sngl" & school != 48,
             index.cond = function(x, y) coef(lm(y ~ x))[1]))


###################################################
### code chunk number 21: Examplot4a
###################################################
print(xyplot(normexam ~ standLRT | school, Exam,
             type = c("g", "p", "r"),
             xlab = "Standardized London Reading Test score",
             ylab = "Normalized exam score",
             subset = sex == "F" & type == "Sngl" & school != 48,
             index.cond = function(x, y) coef(lm(y ~ x))[2]))


###################################################
### code chunk number 22: ExamLmListFS
###################################################
show(ExamFS <- lmList(normexam ~ standLRT | school, Exam,
                      subset = sex == "F" & type == "Sngl" & school != 48))


###################################################
### code chunk number 23: Examplot4cshow
###################################################
plot(confint(ExamFS, pool = TRUE), order = 1)


###################################################
### code chunk number 24: Examplot4c
###################################################
print(plot(confint(ExamFS, pool = TRUE), order = 1))


###################################################
### code chunk number 25: Examplot5
###################################################
print(xyplot(normexam ~ standLRT | school, Exam,
             type = c("g", "p", "r"),
             xlab = "Standardized London Reading Test score",
             ylab = "Normalized exam score",
             subset = sex == "M" & type == "Sngl", layout = c(5,2),
             index.cond = function(x, y) coef(lm(y ~ x))[1]))


###################################################
### code chunk number 26: ExamLmListMS
###################################################
show(ExamMS <- lmList(normexam ~ standLRT | school, Exam,
                      subset = sex == "M" & type == "Sngl"))


###################################################
### code chunk number 27: Examplot5b
###################################################
print(plot(confint(ExamMS, pool = TRUE), order = 1))


###################################################
### code chunk number 28: Examplot6
###################################################
print(xyplot(normexam ~ standLRT | school, Exam, groups = sex,
             type = c("g", "p", "r"),
             xlab = "Standardized London Reading Test score",
             ylab = "Normalized exam score",
             subset = !school %in% c(43, 47) & type == "Mxd",
             index.cond = function(x, y) coef(lm(y ~ x))[1],
             auto.key = list(space = 'top', lines = TRUE,
             columns = 2), layout = c(7,5),
             aspect = 1.2))


###################################################
### code chunk number 29: ExamLmListM
###################################################
show(ExamM <- lmList(normexam ~ standLRT + sex| school, Exam,
                     subset = type == "Mxd" & !school %in% c(43,47,54)))


###################################################
### code chunk number 30: Examplot6b
###################################################
print(plot(confint(ExamM, pool = TRUE), order = 1))


###################################################
### code chunk number 31: Em3
###################################################
(Em3 <- lmer(normexam ~ standLRT + sex + type + (1|school), Exam))


###################################################
### code chunk number 32: Em4
###################################################
(Em4 <- lmer(normexam ~ standLRT + sex + type + (standLRT|school), Exam))


###################################################
### code chunk number 33: EmAnova
###################################################
anova(Em3, Em4)


###################################################
### code chunk number 34: Em5
###################################################
(Em5 <- lmer(normexam ~ standLRT + sex + type + (standLRT + sex|school), Exam))


###################################################
### code chunk number 35: Oxboys
###################################################
str(Oxboys)
system.time(mX1 <- lmer(height ~ age + I(age^2) + I(age^3) + I(age^4) + (age + I(age^2)|Subject),
                       Oxboys))
summary(mX1)
system.time(mX2 <- lmer(height ~ poly(age,4) + (age + I(age^2)|Subject), Oxboys))
summary(mX2)


###################################################
### code chunk number 36: ScotsSec
###################################################
str(ScotsSec)
system.time(mS1 <- lmer(attain ~ sex + (1|primary) + (1|second), ScotsSec))
summary(mS1)


###################################################
### code chunk number 37: sessionInfo
###################################################
toLatex(sessionInfo())


