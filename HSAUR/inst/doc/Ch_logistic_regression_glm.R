### R code from vignette source 'Ch_logistic_regression_glm.Rnw'
### Encoding: ASCII

###################################################
### code chunk number 1: setup
###################################################
rm(list = ls())
s <- search()[-1]
s <- s[-match(c("package:base", "package:stats", "package:graphics", "package:grDevices",
                "package:utils", "package:datasets", "package:methods", "Autoloads"), s)]
if (length(s) > 0) sapply(s, detach, character.only = TRUE)
if (!file.exists("tables")) dir.create("tables")
if (!file.exists("figures")) dir.create("figures")
set.seed(290875)
options(prompt = "R> ", continue = "+  ",
    width = 63, # digits = 4, 
    SweaveHooks = list(leftpar = function() 
        par(mai = par("mai") * c(1, 1.05, 1, 1))))
HSAURpkg <- require("HSAUR")
if (!HSAURpkg) stop("cannot load package ", sQuote("HSAUR"))
rm(HSAURpkg)
a <- Sys.setlocale("LC_ALL", "C")
book <- TRUE
refs <- cbind(c("AItR", "SI", "CI", "ANOVA", "MLR", "GLM", 
                "DE", "RP", "SA", "ALDI", "ALDII", "MA", "PCA", 
                "MDS", "CA"), 1:15)
ch <- function(x, book = TRUE) {
    ch <- refs[which(refs[,1] == x),]
    if (book) {
        return(paste("Chapter~\\\\ref{", ch[1], "}", sep = ""))
    } else {
        return(paste("Chapter~\\\\ref{", ch[2], "}", sep = ""))
    }
}


###################################################
### code chunk number 2: GLM-plasma-plot
###################################################
data("plasma", package = "HSAUR")
layout(matrix(1:2, ncol = 2))
cdplot(ESR ~ fibrinogen, data = plasma)
cdplot(ESR ~ globulin, data = plasma)


###################################################
### code chunk number 3: GLM-plasma-fit1
###################################################
plasma_glm_1 <- glm(ESR ~ fibrinogen, data = plasma, 
                    family = binomial())


###################################################
### code chunk number 4: GLM-plasma-summary-1
###################################################
summary(plasma_glm_1)


###################################################
### code chunk number 5: GLM-plasma-confint
###################################################
ci <- confint(plasma_glm_1)["fibrinogen",]


###################################################
### code chunk number 6: GLM-plasma-confint
###################################################
confint(plasma_glm_1, parm = "fibrinogen")


###################################################
### code chunk number 7: GLM-plasma-confint
###################################################
print(ci) 


###################################################
### code chunk number 8: GLM-plasma-exp
###################################################
exp(coef(plasma_glm_1)["fibrinogen"])


###################################################
### code chunk number 9: GLM-plasma-exp-ci
###################################################
ci <- exp(confint(plasma_glm_1, parm = "fibrinogen"))


###################################################
### code chunk number 10: GLM-plasma-exp-ci
###################################################
exp(confint(plasma_glm_1, parm = "fibrinogen"))


###################################################
### code chunk number 11: GLM-plasma-exp-ci
###################################################
print(ci)


###################################################
### code chunk number 12: GLM-plasma-fit2
###################################################
plasma_glm_2 <- glm(ESR ~ fibrinogen +  globulin, data = plasma, 
                   family = binomial())


###################################################
### code chunk number 13: GLM-plasma-summary-2
###################################################
summary(plasma_glm_2)


###################################################
### code chunk number 14: GLM-plasma-anova-hide
###################################################
plasma_anova <- anova(plasma_glm_1, plasma_glm_2, test = "Chisq")


###################################################
### code chunk number 15: GLM-plasma-anova
###################################################
anova(plasma_glm_1, plasma_glm_2, test = "Chisq")


###################################################
### code chunk number 16: GLM-plasma-predict
###################################################
prob <- predict(plasma_glm_2, type = "response")


###################################################
### code chunk number 17: GLM-plasma-bubble
###################################################
plot(globulin ~ fibrinogen, data = plasma, xlim = c(2, 6), 
     ylim = c(25, 55), pch = ".")
symbols(plasma$fibrinogen, plasma$globulin, circles = prob,
        add = TRUE)


###################################################
### code chunk number 18: GLM-womensrole-fit1
###################################################
data("womensrole", package = "HSAUR")
fm1 <- cbind(agree, disagree) ~ sex + education
womensrole_glm_1 <- glm(fm1, data = womensrole, 
                        family = binomial())


###################################################
### code chunk number 19: GLM-womensrole-summary-1
###################################################
summary(womensrole_glm_1)


###################################################
### code chunk number 20: GLM-womensrole-probfit
###################################################
role.fitted1 <- predict(womensrole_glm_1, type = "response")


###################################################
### code chunk number 21: GLM-plot-setup
###################################################
myplot <- function(role.fitted)  {
    f <- womensrole$sex == "Female"
    plot(womensrole$education, role.fitted, type = "n", 
         ylab = "Probability of agreeing",
         xlab = "Education", ylim = c(0,1))
    lines(womensrole$education[!f], role.fitted[!f], lty = 1)
    lines(womensrole$education[f], role.fitted[f], lty = 2)  
    lgtxt <- c("Fitted (Males)", "Fitted (Females)")
    legend("topright", lgtxt, lty = 1:2, bty = "n")
    y <-  womensrole$agree / (womensrole$agree + 
                              womensrole$disagree)
    text(womensrole$education, y, ifelse(f, "\\VE", "\\MA"), 
         family = "HersheySerif", cex = 1.25)
} 


###################################################
### code chunk number 22: GLM-role-fitted1
###################################################
myplot(role.fitted1)


###################################################
### code chunk number 23: GLM-womensrole-fit2
###################################################
fm2 <- cbind(agree,disagree) ~ sex * education
womensrole_glm_2 <- glm(fm2, data = womensrole, 
                        family = binomial())


###################################################
### code chunk number 24: GLM-womensrole-summary-2
###################################################
summary(womensrole_glm_2)


###################################################
### code chunk number 25: GLM-role-fitted2
###################################################
role.fitted2 <- predict(womensrole_glm_2, type = "response")
myplot(role.fitted2)


###################################################
### code chunk number 26: GLM-role-plot2
###################################################
res <- residuals(womensrole_glm_2, type = "deviance")
plot(predict(womensrole_glm_2), res,
     xlab="Fitted values", ylab = "Residuals", 
     ylim = max(abs(res)) * c(-1,1))
abline(h = 0, lty = 2)


###################################################
### code chunk number 27: GLM-polyps-fit1
###################################################
data("polyps", package = "HSAUR")
polyps_glm_1 <- glm(number ~ treat + age, data = polyps, 
                    family = poisson())


###################################################
### code chunk number 28: GLM-polyps-summary-1
###################################################
summary(polyps_glm_1)


###################################################
### code chunk number 29: GLM-polyp-quasi
###################################################
polyps_glm_2 <- glm(number ~ treat + age, data = polyps, 
                    family = quasipoisson())
summary(polyps_glm_2)


