
library("multcomp")
tol <- sqrt(.Machine$double.eps)
set.seed(29081975)

df <- data.frame(y = rnorm(100), 
                 x = runif(100),  
                 z = runif(100))

### linear model
fam <- gaussian()
lm0 <- glm(y ~ 1, data = df, family = fam)
lm1 <- glm(y ~ x, data = df, family = fam)
lm2 <- glm(y ~ x + z, data = df, family = fam)

gh <- glht(lm2, linfct = c("x = 0", "z = 0"))
stopifnot(abs(anova(lm0, lm2, test = "F")[2, 6] - 
    summary(gh, test = Ftest())$test$pvalue) < tol)
stopifnot(abs(anova(lm0, lm2, test = "Chisq")[2, 5] - 
    summary(gh, test = Chisqtest())$test$pvalue) < tol)

gh <- glht(lm2, linfct = "z = 0")
stopifnot(abs(anova(lm1, lm2, test = "F")[2, 6] - 
    summary(gh, test = Ftest())$test$pvalue) < tol)
stopifnot(abs(anova(lm1, lm2, test = "Chisq")[2, 5] - 
    summary(gh, test = Chisqtest())$test$pvalue) < tol)

### logistic regression
df$y <- factor(df$y < 0)
fam <- binomial()
lm0 <- glm(y ~ 1, data = df, family = fam)
lm1 <- glm(y ~ x, data = df, family = fam)
lm2 <- glm(y ~ x + z, data = df, family = fam)

if (require("lmtest")) {

  gh <- glht(lm2, linfct = c("x = 0", "z = 0"))
  stopifnot(abs(waldtest(lm0, lm2, test = "Chisq")[2, 4] - 
      summary(gh, test = Chisqtest())$test$pvalue) < tol)

  gh <- glht(lm2, linfct = "z = 0")
  stopifnot(abs(waldtest(lm1, lm2, test = "Chisq")[2, 4] -
      summary(gh, test = Chisqtest())$test$pvalue) < tol)
}

