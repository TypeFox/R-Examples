library(VGAM)

# Zelig 4 code:
library(Zelig4)
library(ZeligChoice4)
data(sanction)
sanction$ncost <- factor(sanction$ncost, ordered = TRUE,
                         levels = c("net gain", "little effect",
                                    "modest loss", "major loss"))
z.out <- Zelig4::zelig(ncost ~ mil + coop, model = "ologit", data = sanction)
summary(z.out)
x.out <- Zelig4::setx(z.out, fn = NULL)
set.seed(42)
s.out <- Zelig4::sim(z.out, x = x.out, num = 100)
summary(s.out)

# Zelig 5 code:
data(sanction)
z5 <- zologit$new()
z5
z5$zelig(ncost ~ mil + coop, data = sanction)
z5
z5$setx(coop = 1:3)

z5$setx()

z.out <- z5$zelig.out$z.out[[1]]

set.seed(42)
z5$sim(num = 100)
z5$sim.out
z5$summarize()
z5$cite()

z5 <- zologit$new()
z5
z5$zelig(ncost ~ mil + coop, data = sanction, by = "export")
z5
z5$setx()

z.out <- z5$zelig.out$z.out[[1]]

set.seed(42)
z5$sim(num = 100)
z5$sim.out
z5$summarize()
z5$cite()


fit <- MASS::polr(formula = as.factor(ncost) ~ mil + coop, data = sanction, method = "logistic", Hess = TRUE)
summary(fit)

fit2 <- MASS::polr(formula = ncost ~ mil + coop, data = sanction, method = "logistic", Hess = TRUE)
summary(fit2)
# z5$zelig(list(import ~ coop + cost + target, export ~ coop + cost + target), data = sanction)

