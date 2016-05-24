library(gnm)
set.seed(1)
DNase1 <- subset(DNase, Run == 1)

fm3DNase1.2 <- gnm(density ~ -1 +
                   Mult(1, Inv(Const(1) + Exp(1 + Mult(offset(-log(conc)),
                                                       Inv(1))))),
                   start = c(NA, 0, 1), data = DNase1, trace = TRUE)
coef(fm3DNase1.2)
