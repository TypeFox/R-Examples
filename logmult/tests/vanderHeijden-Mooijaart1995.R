## van der Heijden & Mooijaart (1995), Table 2c, p. 23
# See also ?hmskew

library(logmult)
data(ocg1973)

# 5:1 is here to take "Farmers" as reference category (angle 0)
model <- hmskew(ocg1973[5:1, 5:1], weighting="uniform", start=NA)
model
ass <- model$assoc.hmskew

# First column of the table
round(ass$row[,,1] * sqrt(ass$phi[1,1]), d=2)[5:1,]

summary(model)

# First score for Farmers is slightly different from the original article
stopifnot(isTRUE(all.equal(round(ass$row[,,1] * sqrt(ass$phi[1,1]), d=2)[5:1,],
                                 matrix(c(-0.08, -0.2, -0.23, -0.11, 0.61,
                                           0.34,  0.3, -0.13, -0.51, 0), 5, 2),
                          check.attributes=FALSE)))

# Right part of the table
round(ass$phi[1] * (ass$row[,2,1] %o% ass$row[,1,1] - ass$row[,1,1] %o% ass$row[,2,1]), d=3)[5:1, 5:1]

# Plot
plot(model, coords="cartesian")

# Test anova
indep <- gnm(Freq ~ O + D, data=ocg1973, family=poisson)
symm <- gnm(Freq ~ O + D + Symm(O, D), data=ocg1973, family=poisson)
anova(indep, symm, model, test="LR")
