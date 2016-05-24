
library("mlt")
library("lattice")

data("faithful")

aic <- numeric(20)

for (o in 2:(length(aic) + 1)) {
Bs <- Bernstein_basis(numeric_var("waiting", support = range(faithful$waiting) + c(-5, 5)),
                      order = o, ui = "incre")
m <- ctm(Bs)
mod <- mlt(m, data = faithful)
yp <- mkgrid(mod, 50)[["waiting"]]

aic[o - 1] <- AIC(mod)

pd <- data.frame(waiting = yp)
pd$p <- predict(mod, q = yp, type = "distribution")

plot(p ~ waiting, data = pd,
     col = "red", pch = 21, main = paste("order", o, "aic", aic[o - 1]))
lines(ecdf(faithful$waiting))

}

plot(aic)

o <- which.min(aic) + 1
Bs <- Bernstein_basis(numeric_var("waiting", support = range(faithful$waiting) + c(-5, 5)),
                      order = o, ui = "incre")
m <- ctm(Bs)
mod <- mlt(m, data = faithful)

abline(h = AIC(mod))

pd$d <- predict(mod, q = yp, type = "density")

plot(d ~ waiting, data = pd, type = "l", col = "red", lwd = 3)
lines(density(faithful$waiting))
lines(rug(faithful$waiting))
abline(h = 0)

p <- 1:99 / 100
q <- predict(mod, p = p, K = 100, type = "quantile")

plot(p, q)
lines(p, quantile(faithful$waiting, p))

Bs <- Bernstein_basis(numeric_var("waiting", support = range(faithful$waiting) + c(-5, 5)),
                      order = o, ui = "incre")
m <- ctm(Bs)
mod <- mlt(m, data = faithful)

# H1 <- mod$optim(coef(mod), hessian = TRUE)$hessian
H2 <- mod$hessian(coef(mod), weights(mod))

X <- model.matrix(m, faithful)
Xprime <- model.matrix(m, faithful, deriv = c(waiting = 1))
w <- drop((Xprime %*% coef(mod))^2)
H3 <- crossprod(X) + crossprod(Xprime * w, Xprime)
max(abs(H3 - H2))

cov2cor(vcov(mod))

if (FALSE) {
library("multcomp") ### since 1.0-3

mp <- parm(coef(mod), vcov(mod))
y <- mkgrid(mod, 30)$waiting
g <- glht(mp, linfct = model.matrix(mod$model,
    data = data.frame(waiting = y)))

mc <- confint(g)
umc <- confint(g, calpha = qnorm(.975))
p <- mod$model$todistr$p
plot(y, p(mc$confint[, "Estimate"]), type = "l")
lines(y, p(mc$confint[, "lwr"]))
lines(y, p(mc$confint[, "upr"]))
lines(y, p(umc$confint[, "lwr"]))
lines(y, p(umc$confint[, "upr"]))

library("survival")
cm <- coxph(Surv(waiting, rep(TRUE, nrow(faithful))) ~ 1, data = faithful)
plot(survfit(cm))
lines(y, 1 - p(mc$confint[, "Estimate"]), col = "red")
lines(y, 1 - p(mc$confint[, "lwr"]), col = "red")
lines(y, 1 - p(mc$confint[, "upr"]), col = "red")
}
