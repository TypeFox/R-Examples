data("data_122")

## Model formula:
fr1 <- y ~ x1 + x2 | x3 + f1 | 0 | x4 + x5 | z1 + z2 | z3 + f2
##     ^   ^^^^^^^   ^^^^^^^   ^   ^^^^^^^   ^^^^^^^   ^^^^^^^
##     y     u11       u12    u13    u14       u22       u24

m1 <- egame122(fr1, data = data_122)
summary(m1)

## Dummy specification of the dependent variable
fr2 <- update(Formula(fr1), a1 + a2 ~ .)
m2 <- egame122(fr2, data = data_122)
summary(m2)
