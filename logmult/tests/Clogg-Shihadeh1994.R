## Clogg & Shihadeh (1994), Tables 3.5a and b (p. 55-61)
# See also ?rc

library(logmult)
data(gss88)
model <- rc(gss88, start=NA)

# Unweighted scores
summary(model, weighting="none")
# Marginally weighted scores
summary(model, weighting="marginal")
# Uniformly weighted scores
summary(model, weighting="uniform")

unweighted <- assoc(model, weighting="none")
marginal <- assoc(model, weighting="marginal")
uniform <- assoc(model, weighting="uniform")

stopifnot(all.equal(round(c(unweighted$row, unweighted$row * sqrt(unweighted$phi[1]),
                            marginal$row, marginal$row * sqrt(marginal$phi[1]),
                            uniform$row, uniform$row * sqrt(uniform$phi[1])), 2),
                    c(0.57, 0.56, 0.15, -0.01, -0.23, 0.04, -0.42, -0.2, -0.19, -0.16, -0.1,
                      1.61, 1.6, 0.42, -0.04, -0.66, 0.11, -1.19, -0.56, -0.54, -0.46, -0.29,
                      1.78, 1.78, 0.37, -0.17, -0.91, 0, -1.53, -0.78, -0.76, -0.66, -0.47,
                      1.53, 1.52, 0.32, -0.14, -0.78, 0, -1.31, -0.67, -0.65, -0.57, -0.4,
                      1.88, 1.87, 0.48, -0.04,  -0.77, 0.12, -1.38, -0.65, -0.63, -0.53, -0.34, 
                      1.87, 1.87, 0.48, -0.04, -0.77, 0.12, -1.38, -0.65, -0.63, -0.53, -0.34)))
stopifnot(all.equal(round(c(unweighted$col, unweighted$col * sqrt(unweighted$phi[1]),
                            marginal$col, marginal$col * sqrt(marginal$phi[1]),
                            uniform$col, uniform$col * sqrt(uniform$phi[1])), 2),
                    c(-0.56, -0.41, -0.08, 0.07, 0.41, 0.58,
                      -1.6, -1.16, -0.24, 0.19, 1.16, 1.65,
                      -1.89, -1.39, -0.33, 0.16, 1.29, 1.85,
                      -1.62, -1.19, -0.28, 0.14, 1.1, 1.59,
                      -1.38, -1, -0.21, 0.16, 1, 1.43,
                      -1.37, -1, -0.21, 0.16, 1, 1.42)))

# Test anova
indep <- gnm(Freq ~ Occupation + Years.of.Schooling, data=gss88, family=poisson)
anova(indep, unweighted, test="LR")
