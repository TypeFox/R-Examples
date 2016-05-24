data("data_123")

## Model formula:
f1 <- y ~ x1 + x2 | 0 | x3 | x4 + x5 | 0 | x6 | x7 | x8
##    ^   ^^^^^^^   ^   ^^   ^^^^^^^   ^   ^^   ^^   ^^
##    y     u11    u13  u15    u16    u23  u25  u26  u36

m1 <- egame123(f1, data = data_123, link = "probit", type = "private")
summary(m1)

## Dummy specification of the dependent variable
f2 <- update(Formula(f1), a1 + a2 + a3 ~ .)
m2 <- egame123(f2, data = data_123, link = "probit", type = "private")
summary(m2)
