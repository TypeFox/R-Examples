
library("mlt")
library("survival")

data("GBSG2", package = "TH.data")

xvar <- names(GBSG2)
xvar <- xvar[!(xvar %in% c("time", "cens"))]
GBSG2$y <- with(GBSG2, Surv(time, cens))

fm <- as.formula(paste("Surv(time, cens) ~ ", paste(xvar, collapse = "+")))
cmod <- coxph(fm, data = GBSG2)

order <- 10
by <- Bernstein_basis(numeric_var("y", support = c(0, max(GBSG2$time))), order = order,
                      ui = "incre")
bx <- as.basis(as.formula(paste("~", paste(xvar, collapse = "+"))), data = GBSG2,
               remove_intercept = TRUE)

m <- ctm(by, shift = bx, todist = "MinEx")

mod <- mlt(m, data = GBSG2, scale = TRUE, check = FALSE)

n <- names(coef(cmod))
cf <- coef(mod)[n]
v <- vcov(mod)[n, n]
coef(cmod) / cf
diag(vcov(cmod)) / diag(v)
range(vcov(cmod) / v)
