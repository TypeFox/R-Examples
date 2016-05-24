library(VGAM)

## Results don't match: Zelig 4 seems to be using the logit inverse link in the probit model

# Zelig 4 code:
library(Zelig4)
library(ZeligChoice4)
data(sanction)
sanction$ncost <- factor(sanction$ncost, ordered = TRUE,
                         levels = c("net gain", "little effect",
                                    "modest loss", "major loss"))
z.out <- Zelig4::zelig(ncost ~ mil + coop, model = "oprobit", data = sanction)
summary(z.out)
x.out <- Zelig4::setx(z.out, fn = NULL)
set.seed(42)
s.out <- Zelig4::sim(z.out, x = x.out, num = 5)
summary(s.out)

# Zelig 5 code:
data(sanction)
z5 <- zoprobit$new()
z5
z5$zelig(ncost ~ mil + coop, data = sanction)
z5
z5$setrange(sanction = 1)
z5
z5$sim(num = 100)
z5
z5$setx()

set.seed(42)
z5$sim(num = 5)
z5$sim.out
z5$summarize()
z5$cite()
