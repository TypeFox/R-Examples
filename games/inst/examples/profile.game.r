data("student_offers")

## A model that does not converge to global max
f1 <- offer + accept ~ gender1 | gender2
m1 <- ultimatum(f1, maxOffer = 100, data = student_offers, s2 = 1)

p1 <- profile(m1)  ## Issues warning
plot(p1)

## Refit model with better starting values
m2 <- ultimatum(f1, maxOffer = 100, data = student_offers, s2 = 1, profile = p1)
p2 <- profile(m2)
plot(p2)

logLik(m1)
logLik(m2)  ## Improved
