data("war1800")

## Model formula:
f1 <- esc + war ~ s_wt_re1 + revis1 | 0 | regime1 | balanc + regime2
##    ^^^^^^^^^   ^^^^^^^^^^^^^^^^^   ^   ^^^^^^^   ^^^^^^^^^^^^^^^^
##        y              u11         u13    u14           u24

m1 <- egame12(f1, data = war1800)
summary(m1)

m2 <- egame12(f1, data = war1800, link = "logit")
summary(m2)

m3 <- egame12(f1, data = war1800, subset = year >= 1850)
summary(m3)

m4 <- egame12(f1, data = war1800, boot = 10)
summary(m4)
summary(m4, useboot = FALSE)

## Estimating scale parameters under fixed utilities
utils <- c(-1, 0, -1.4, 0.1)
m5 <- egame12(esc + war ~ 1, data = war1800, fixedUtils = utils)
summary(m5)

m6 <- egame12(esc + war ~ 1, data = war1800, fixedUtils = utils, sdByPlayer = TRUE)
summary(m6)

## Estimating scale parameters with regressors
m7 <- egame12(f1, data = war1800, sdformula = ~ balanc - 1)
summary(m7)

## Using a factor outcome
y <- ifelse(war1800$esc == 1, ifelse(war1800$war == 1, "war", "cap"), "sq")
war1800$y <- factor(y, levels = c("sq", "cap", "war"))
f2 <- update(Formula(f1), y ~ .)

m8 <- egame12(f2, data = war1800)
summary(m8)
