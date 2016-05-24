
library("multcomp")
set.seed(290875)

testdata <- data.frame(y = rnorm(21), 
                       f1 <- factor(c(rep(c("A", "B", "C"), 7))),
                       f2 <- factor(c(rep("D", 10), rep("E", 11))),
                       x <- rnorm(21))

# one-way ANOVA
coef(amod <- aov(y ~ f1, data = testdata))
glht(amod, linfct = mcp(f1 = "Dunnett"))

# and a continuous covariable: ANCOVA
coef(lmod <- lm(y ~ f1 + x, data = testdata))
glht(lmod, linfct = mcp(f1 = "Dunnett"))

# ANCOVA with an additional factor as covariable
coef(lmod <- lm(y ~ f1 + f2 + x, data = testdata))
glht(lmod, linfct = mcp(f1 = "Dunnett"))

# and with interaction terms
coef(lmod <- lm(y ~ f1 + f2 + f2:f1 + x, data = testdata))
glht(lmod, linfct = mcp(f1 = "Dunnett"))

# with contrasts as expressions
glht(lmod, linfct = mcp(f1 = c("B - A = 0", "C - A = 0")))

tmp <- multcomp:::chrlinfct2matrix(c(l1 = "x1 - x2 = 2", 
                                      l2 = "x2 + 3 * x3 = 1"), 
                                      paste("x", 1:3, sep = ""))

stopifnot(max(abs(tmp$K - rbind(c(1, -1, 0), c(0, 1, 3)))) < sqrt(.Machine$double.eps))
stopifnot(max(abs(tmp$m - c(2, 1))) < sqrt(.Machine$double.eps))

### coef.survreg and vcov.survreg need special tuning
### thx to Z for pointing this out
if (require("survival")) {
    smod <- survreg(Surv(futime, fustat) ~ ecog.ps + rx, 
                    data = ovarian, dist = 'weibull')
    K <- diag(length(coef(smod)))
    rownames(K) <- names(coef(smod))
    glht(smod, linfct = K)
}

### new `means' comparisons
amod <- aov(weight ~ dose + gesttime + number, data = litter)
confint(glht(amod, linfct = mcp(dose = "Means")))
