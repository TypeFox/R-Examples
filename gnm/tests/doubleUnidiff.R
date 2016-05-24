library(gnm)
set.seed(1)

doubleUnidiff <- gnm(Freq ~ election*vote + election*class*religion +
                     Mult(Exp(election), religion:vote) +
                     Mult(Exp(election), class:vote),
                     family = poisson, data = cautres)

print(doubleUnidiff$deviance, digits=10)
print(doubleUnidiff$df)
