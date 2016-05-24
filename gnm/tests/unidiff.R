library(gnm)
set.seed(1)

unidiff <- gnm(Freq ~ educ*orig + educ*dest +
               Mult(Exp(educ), orig:dest), family = poisson,
               data = yaish, subset = (dest != 7))

print(unidiff$deviance, digits = 10)
print(unidiff$df)

getContrasts(unidiff, grep("[.]educ", names(coef(unidiff))))
