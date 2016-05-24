## Examples

library(bootStepAIC)

## lm() Example ##
n <- 100
x1 <- runif(n, -4, 4)
x2 <- runif(n, -4, 4)
x3 <- runif(n, -4, 4)
x4 <- runif(n, -4, 4)
x5 <- runif(n, -4, 4)
x6 <- runif(n, -4, 4)
x7 <- factor(sample(letters[1:3], n, rep = TRUE))
y <- 5 + 3 * x1 + 2 * x2 - 1.5 * x3 - 0.8 * x4 + rnorm(n, sd = 2.5)
data <- data.frame(y, x1, x2, x3, x4, x5, x6, x7)
rm(n, x1, x2, x3, x4, x5, x6, x7, y)

lmFit <- lm(y ~ ., data = data)
boot.stepAIC(lmFit, data)


## glm() Example ##
n <- 200
x1 <- runif(n, -3, 3)
x2 <- runif(n, -3, 3)
x3 <- runif(n, -3, 3)
x4 <- runif(n, -3, 3)
x5 <- factor(sample(letters[1:2], n, rep = TRUE))
eta <- 0.1 + 1.6 * x1 - 2.5 * as.numeric(as.character(x5) == levels(x5)[1])
y1 <- rbinom(n, 1, plogis(eta))
y2 <- rbinom(n, 1, 0.6)
data <- data.frame(y1, y2, x1, x2, x3, x4, x5)
rm(n, x1, x2, x3, x4, x5, eta, y1, y2)

glmFit1 <- glm(y1 ~ x1 + x2 + x3 + x4 + x5, family = binomial, data = data)
glmFit2 <- glm(y2 ~ x1 + x2 + x3 + x4 + x5, family = binomial, data = data)

boot.stepAIC(glmFit1, data, B = 50)
boot.stepAIC(glmFit2, data, B = 50)


## aov() Example ##
quine.hi <- aov(log(Days + 2.5) ~ .^4, quine)
quine.nxt <- update(quine.hi, . ~ . - Eth:Sex:Age:Lrn)
boot.stepAIC(quine.nxt, data = quine, B = 100,
    scope = list(upper = ~ Eth * Sex * Age * Lrn, lower = ~ 1))


## glm.nb() Example ##
quine.nb <- glm.nb(Days ~ .^4, data = quine)
boot.stepAIC(quine.nb, data = quine)


## polr() Example ##
house <- housing[rep(1:nrow(housing), housing$Freq), ]
house.plr <- polr(Sat ~ Infl + Type + Cont, data = house)
boot.stepAIC(house.plr, data = house)


## survreg() Example ##
require(survival)
survFit <- survreg(Surv(futime, fustat) ~ ecog.ps + rx, ovarian, dist = "weibull")
boot.stepAIC(survFit, data = ovarian)


## coxph() Example ##
coxFit <- coxph(Surv(futime, fustat) ~ ecog.ps + rx, ovarian)
boot.stepAIC(coxFit, data = ovarian, B = 20)
