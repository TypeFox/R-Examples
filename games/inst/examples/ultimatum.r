data(data_ult)

## Model formula:
f1 <- offer + accept ~ x1 + x2 + x3 + x4 + w1 + w2 | z1 + z2 + z3 + z4 + w1 + w2
##                     ^^^^^^^^^^^^^^^^^^^^^^^^^^^   ^^^^^^^^^^^^^^^^^^^^^^^^^^^
##                                  R1                              R2

m1 <- ultimatum(f1, data = data_ult, maxOffer = 15)
summary(m1)

## Estimating offer size only
f2 <- update(Formula(f1), offer ~ .)
m2 <- ultimatum(f2, data = data_ult, maxOffer = 15, outcome = "offer")
summary(m2)

## Fixing scale terms
m3 <- ultimatum(f1, data = data_ult, maxOffer = 15, s1 = 5, s2 = 1)
summary(m3)
