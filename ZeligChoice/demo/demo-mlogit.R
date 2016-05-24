library(VGAM)

# Zelig 4 code:
library(Zelig4)
library(ZeligChoice4)
data(mexico)

z.out1 <- Zelig4::zelig(as.factor(vote88) ~ pristr + othcok + othsocok,
                        model = "mlogit", data = mexico)
summary(z.out1)
x.weak <- Zelig4::setx(z.out1, pristr = 1)
x.strong <- Zelig4::setx(z.out1, pristr = 3)
x.out <- Zelig4::setx(z.out1)
set.seed(42)
s.out1 <- Zelig4::sim(z.out1, x = x.out)
summary(s.out1)

v <- VGAM::vglm(formula = as.factor(vote88) ~ pristr + othcok + 
                othsocok, data = mexico, family = "multinomial")

# Zelig 5 code:
data(mexico)
z5 <- zmlogit$new()
z5
z5$zelig(as.factor(vote88) ~ pristr + othcok + othsocok, data = mexico)
z5
z5$setx()
set.seed(42)
z5$sim(num = 1000)
z5$sim.out
z5$summarize()
z5$cite()

# z5$zelig(list(import ~ coop + cost + target, export ~ coop + cost + target), data = sanction)

