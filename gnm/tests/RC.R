library(gnm)
set.seed(1)

mentalHealth$MHS <- C(mentalHealth$MHS, treatment)
mentalHealth$SES <- C(mentalHealth$SES, treatment)
RC1model <- gnm(count ~ SES + MHS +
                Mult(-1 + SES, -1 + MHS),
                family = poisson, data = mentalHealth)

print(RC1model$deviance, digits = 10)
print(RC1model$df)

rowProbs <- with(mentalHealth, tapply(count, SES, sum) / sum(count))
colProbs <- with(mentalHealth, tapply(count, MHS, sum) / sum(count))
mu <- getContrasts(RC1model, pickCoef(RC1model, "[.]SES"),
                   ref = rowProbs, scaleRef = rowProbs,
                   scaleWeights = rowProbs)
nu <- getContrasts(RC1model, pickCoef(RC1model, "[.]MHS"),
                   ref = colProbs, scaleRef = colProbs,
                   scaleWeights = colProbs)
all.equal(sum(mu$qv[,1] * rowProbs), 0)
all.equal(sum(nu$qv[,1] * colProbs), 0)
all.equal(sum(mu$qv[,1]^2 * rowProbs), 1)
all.equal(sum(nu$qv[,1]^2 * colProbs), 1)
