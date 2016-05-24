# some ingredients
x.bar <- mean(~ BodyTemp, data = BodyTemp50); x.bar
sd <- sd (~ BodyTemp, data = BodyTemp50); sd
n <- nrow (BodyTemp50); n
se <- sd/sqrt(n); se
# test statistic
t <- (x.bar - 98.6) / se; t
# 2-sided p-value
2 * pt(- abs(t), df = 49)

